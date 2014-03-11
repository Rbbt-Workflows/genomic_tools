require 'rbbt/entity/genomic_mutation'

module GenomicTools::ExomeSequencing
  extend Workflow

  helper :organism do
    "Hsa/jan2013"
  end

  helper :watson do
    true
  end

  helper :principal_isoforms do
    Appris::PRINCIPAL_ISOFORMS
  end

  task :MISSING_SAMPLE => nil do
    raise "Sample #{ name } is missing"
  end

  dep :MISSING_SAMPLE
  task :BED => :text do
  end

  desc "This is not implemented yet. Should use Miriams RuBioSeq"
  dep :BED
  extension :vcf
  task :VCF => :text do
    raise "Cannot produce VCF from BED just yet"
  end

  desc "Reads the VCF file and filters the variants according to the specified criteria. Currently only a quality threshold."
  dep :VCF 
  input :threshold, :integer, "Quality threshold", 200
  task :GenomicMutations => :array do |threshold|
    vcf_file = step(:VCF).path
    stream = GenomicMutation::VCF.open_stream(vcf_file.open)

    mutations = []
    while line = stream.gets
      next if line[0] == "#"
      mutation, _id, quality = line.split("\t")
      mutations << mutation if quality.to_f > threshold
    end
    mutations
  end

  desc "Selects variants with 'coding relevance' which either produce a\
  non-synonymous variation in an isoform, or that may affect splicing (even if synonymous)"
  dep :GenomicMutations 
  extension :list
  task :CodingRelevantMutations => :array do 
    mutations = step(:GenomicMutations).load
    GenomicMutation.setup(mutations, name.to_s, organism, watson)

    mutations.select_by(:relevant?)
  end

  desc "Find variants that potentially affect splicing of some isoform"
  dep :GenomicMutations
  extension :list
  task :SplicingMutations => :array do 
    genomic_mutation = step(:GenomicMutations).load
    GenomicMutation.setup(genomic_mutation, name.to_s, organism, watson)

    genomic_mutation.select_by(:transcripts_with_affected_splicing){|t| t and t.any?}
  end

  desc "Find transcripts potentially affected by splicing mutations"
  dep :SplicingMutations
  extension :list
  task :TranscriptsWithAffectedSplicing=> :tsv do 
    genomic_mutation = step(:SplicingMutations).load
    GenomicMutation.setup(genomic_mutation, name.to_s, organism, watson)

    tsv = Misc.process_to_hash(genomic_mutation){|muts| muts.transcripts_with_affected_splicing }
    TSV.setup(tsv, :key_field => "Genomic Mutation", :fields => ["Ensembl Transcript ID"], :type => :flat, :namespace => organism)
  end

  desc "Find potentially relevant variants"
  dep :CodingRelevantMutations
  dep :SplicingMutations
  extension :list
  task :RelevantMutations => :array do 
    (step(:CodingRelevantMutations).load + 
     step(:SplicingMutations).load).uniq
  end

  desc "Associates variants to their effect (amino acid substitutions) on their different isoforms"
  dep :CodingRelevantMutations 
  extension :tsv
  task :MutatedIsoforms => :tsv do 
    genomic_mutation = step(:CodingRelevantMutations).load
    GenomicMutation.setup(genomic_mutation, name.to_s, organism, watson)

    tsv = Misc.process_to_hash(genomic_mutation){|muts| muts.mutated_isoforms }

    TSV.setup(tsv, :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat, :namespace => organism)

    tsv
  end

  desc "Protein domains (InterPro) truncated or mutated"
  dep :MutatedIsoforms 
  extension :tsv
  task :AffectedDomains => :tsv do 
    all_mis = step(:MutatedIsoforms).path.tsv(:unnamed => true).values.compact.flatten
    tsv = TSV.setup(all_mis, :key_field => "Mutated Isoform", :fields => [], :type => :double, :namespace => organism)

    tsv.add_field "InterPro Domain" do |mi, values|
      mi.affected_domains
    end

    tsv.add_field "InterPro Domain Name" do |mi, values|
      domains = values["InterPro Domain"]
      domains.any? ? domains.name : []
    end

    tsv
  end

  desc "Reports mutated isoforms only on principal isoforms according to Appris"
  dep :MutatedIsoforms
  extension :tsv
  task :MutatedPrincipalIsoforms => :tsv do 
    tsv = step(:MutatedIsoforms).load
    tsv.process "Mutated Isoform" do |mis|
      mis.select{|mi| principal_isoforms.include? mi.protein}
    end
  end

  desc "List of genes with affected principal isoforms"
  dep :MutatedPrincipalIsoforms
  task :AffectedGenes => :array do 
    mis = MutatedIsoform.setup(step(:MutatedPrincipalIsoforms).load.values.compact.flatten.uniq, organism)
    mis.select_by(:non_synonymous).protein.compact.uniq.gene.compact.uniq
  end

  desc "Filters the mutated isoforms to retain only those considered damaged (details to come)"
  dep :MutatedPrincipalIsoforms
  extension :tsv
  task :DamagedPrincipalIsoforms => :tsv do 
    tsv = step(:MutatedPrincipalIsoforms).load
    all_mis = MutatedIsoform.setup(tsv.values.compact.flatten.uniq, organism)

    damaged_mis = all_mis.select_by(:damaged?)
    tsv.process "Mutated Isoform" do |mis|
      mis & damaged_mis
    end
  end

  dep :MutatedIsoforms
  task :AffectedProteinFeaturesUniProt => :tsv do
    mis = step(:Mutatedso).load
    Workflow.require_workflow "Structure"

    miss_sense = mis.select_by(:consequence, "MISS-SENSE")
    FileUtils.cp Structure.job(:mutated_isoform_annotation, name, :mutated_isoforms => miss_sense).run(true).file(:uniprot), path

    nil
  end

  desc "Variant Consequence"
  dep :GenomicMutations
  dep :SplicingMutations
  dep :MutatedIsoforms
  dep :MutatedPrincipalIsoforms
  dep :TranscriptsWithAffectedSplicing
  extension :tsv
  task :VariantConsequence => :tsv do 
    variants = step(:GenomicMutations).load
    variant_splicing_transcripts = step(:TranscriptsWithAffectedSplicing).load
    splicing = Set.new(step(:SplicingMutations).load)
    mis = step(:MutatedIsoforms).load
    mpis = step(:MutatedPrincipalIsoforms).load
    tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Consequence", "Ensembl Gene ID"], :type => :list)

    enst2ensg = Organism.gene_transcripts(organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Gene ID"], :type => :single, :unnamed => true, :persist => true

    variants.each do |variant|
      consecuence, gene = case
                    when (mpis.include?(variant) and (_mpis = mpis[variant]) and _mpis.any?)
                      [_mpis.first.consequence, _mpis.first.protein.gene]
                    when (mis.include?(variant) and (_mis = mpis[variant]) and _mis.any?)
                      [_mis.first.consequence, _mis.first.protein.gene]
                    else
                      "SILENT"
                    end
      if %w(SYNONYMOUS SILENT).include? consecuence and splicing.include? variant
        consequence = "Exon Junction"
        gene = enst2ensg[variant_splicing_transcripts[variant].first]
      end

      tsv[variant] = [consecuence, gene]
    end

    tsv
  end

  dep :MutatedPrincipalIsoforms
  task :AffectedProteinFeatures => :tsv do 
    Workflow.require_workflow "Structure"
    mpis = step(:MutatedPrincipalIsoforms).load.values.compact.flatten
    job = Structure.job(:mutated_isoform_annotation, name, :mutated_isoforms => mpis)
    job.run

    tsv = nil
    job.files.each do |file|
      log "Attaching #{ file } to #{Misc.fingerprint tsv} from #{job.file(file).find}"
      Open.write(self.file(file), job.file(file).open)
      if tsv.nil?
        tsv = job.file(file).tsv
      else
        tsv.attach job.file(file).tsv
      end
    end

    tsv
  end

  dep :MutatedPrincipalIsoforms
  task :CloseProteinFeatures => :tsv do 
    Workflow.require_workflow "Structure"
    mpis = step(:MutatedPrincipalIsoforms).load.values.compact.flatten
    job = Structure.job(:mutated_isoform_neighbour_annotation, name, :mutated_isoforms => mpis)
    job.run

    tsv = nil
    job.files.each do |file|
      log "Attaching #{ file } to #{Misc.fingerprint tsv} from #{job.file(file).find}"
      Open.write(self.file(file), job.file(file).open)
      if tsv.nil?
        tsv = job.file(file).tsv
      else
        tsv.attach job.file(file).tsv
      end
    end

    tsv
  end

  dep :MutatedPrincipalIsoforms
  task :AffectedComplexInterfaces => :tsv do 
    mpis = step(:MutatedPrincipalIsoforms).load.values.compact.flatten
    job = Structure.job(:mutated_isoform_interface_residues, name, :mutated_isoforms => mpis)
    job.run
  end

  task :ChannelCounts => :tsv do
    raise "NOT IMPLEMENTED"
  end

  dep :MutatedPrincipalIsoforms, :AffectedDomains, :AffectedGenes, :RelevantMutations, :VariantConsequence, :CloseProteinFeatures, :AffectedProteinFeatures, :AffectedComplexInterfaces
  task :all => :string do
    "DONE"
  end
end

