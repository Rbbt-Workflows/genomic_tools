require 'rbbt/entity/genomic_mutation'
module GenomicTools::ExomeSequencing
  extend Workflow

  helper :organism do
    "Hsa/jan2013"
  end

  helper :watson do
    true
  end

  helper :appris_principal_isoforms do
    Appris::PRINCIPAL_ISOFORMS
  end

  task :MISSING_SAMPLE => nil do
    raise "Sample #{ name } is missing"
  end

  dep :MISSING_SAMPLE
  task :BED => :text do
  end

  desc "This is not implemented yet. Should use Miriams RuBioSeq"
  #dep :BED
  extension :vcf
  task :VCF => :text do
    raise "Cannot produce VCF from BED just yet"
  end

  desc "Reads the VCF file and filters the variants according to the specified criteria. Currently only a quality threshold."
  dep :VCF 
  input :threshold, :integer, "Quality threshold", 200
  task :mutations => :array do |threshold|
    step(:VCF).join
    mis = Sequence.job(:genomic_mutations, name, :vcf_file => step(:VCF).path, :threshold => threshold).run 
    mis
  end

  dep :mutations
  task :isoforms => :array do 
    job = Sequence.job(:mutated_isoforms, name, :mutations => step(:mutations).load)
    job.run
    all_mis = []
    TSV.traverse job  do |k,mis|
      next if mis.nil? or mis.empty?
      all_mis.concat mis.flatten 
    end
    all_mis
  end

  dep :isoforms
  task :principal_isoforms => :array do 
    principal_isoforms = self.appris_principal_isoforms
    s = TSV.traverse step(:isoforms), :type => :array, :into => :stream do |mi|
      protein, _sep, change = mi.partition ":"
      next unless principal_isoforms.include? protein
      mi
    end
    s.read.split("\n")
  end

  dep :principal_isoforms
  task :affected_genes => :array do 
    protein_gene = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"], 
      :unnamed => true, :persist => true
    TSV.traverse step(:principal_isoforms), :into => :stream, :type => :array do |mi|
      protein, _sep, change = mi.partition ":"
      protein_gene[protein]
    end
  end

  annotation_tasks = []
  Structure.task_info(:annotated_variants)[:input_options][:database][:select_options].each do |database|
    task_name = [database,"annotations"] * "_"
    task_name = task_name.to_sym

    dep do |jobname,options|
      principal_isoforms_job = GenomicTools::ExomeSequencing.job(:principal_isoforms, jobname, options)
      if not principal_isoforms_job.done?
        principal_isoforms_job.run 
      end
      options = options.merge :mutated_isoforms => principal_isoforms_job, :database => database
      Structure.job(:annotated_variants, jobname, options)
    end
    task task_name => :tsv do 
      step(:annotated_variants).load
    end
    annotation_tasks << task_name

    task_name = [database,"neighbour_annotations"] * "_"
    task_name = task_name.to_sym

    dep do |jobname,options|
      principal_isoforms_job = GenomicTools::ExomeSequencing.job(:principal_isoforms, jobname, options)
      if not principal_isoforms_job.done?
        principal_isoforms_job.run 
      end
      options = options.merge :mutated_isoforms => principal_isoforms_job, :database => database
      Structure.job(:annotated_variant_neighbours, jobname, options)
    end
    task task_name => :tsv do 
      step(:annotated_variant_neighbours).load
    end
    annotation_tasks << task_name
  end
  
  dep do |jobname,options|
    job = GenomicTools::ExomeSequencing.job(:principal_isoforms, jobname, options)
    principal_isoforms = job.run
    job.join

    options = options.merge :mutated_isoforms => principal_isoforms
    Structure.job(:variant_interfaces, jobname, options)
  end
  task :interfaces => :tsv do
    step(:variant_interfaces).load
  end

  dep :interfaces
  dep :affected_genes
  annotation_tasks.each do |dep_name|
    dep dep_name
  end
  task :all => :string do
    "DONE"
  end

  dep :VCF
  dep :isoforms
  task :enhanced_VCF => :tsv do
    mutated_isoforms = step(:isoforms).step(:mutated_isoforms).load

    TSV.traverse step(:VCF), :type => :array, :into => :stream do |line|
      next if line =~/^#/
      chr, pos, _id, ref, alt = line.split("\t")
      chr.sub!(/^chr/,'')

      position, alt = Misc.correct_vcf_mutation(pos.to_i, ref, alt)
      mutation = [chr, position, alt ] * ":"

      isoforms = mutated_isoforms[mutation]
      isoforms = [] if isoforms.nil?
      line << "\t" << isoforms * "|"
    end
  end
end
