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
  dep :BED
  extension :vcf
  task :VCF => :text do
    raise "Cannot produce VCF from BED just yet"
  end

  desc "Reads the VCF file and filters the variants according to the specified criteria. Currently only a quality threshold."
  dep :VCF 
  input :threshold, :integer, "Quality threshold", 200
  task :genomic_mutations => :array do |threshold|
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

  dep do |jobname,options|
    genomic_mutations_job = GenomicTools::ExomeSequencing.job(:genomic_mutations, jobname, options)
    genomic_mutations_job.run

    job = Sequence.job(:mutated_isoforms, jobname, 
                 options.merge(:organism => genomic_mutations_job.organism, :mutations => genomic_mutations_job))
  end
  task :isoforms => :array do 
    TSV.traverse step(:mutated_isoforms), :into => :stream  do |k,mis|
      mis * "\n"
    end
  end

  dep :isoforms
  task :principal_isoforms => :array do 
    principal_isoforms = self.appris_principal_isoforms
    TSV.traverse step(:isoforms), :type => :array, :into => :stream do |mi|
      protein, _sep, change = mi.partition ":"
      next unless principal_isoforms.include? protein
      mi
    end
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
      principal_isoforms = GenomicTools::ExomeSequencing.job(:principal_isoforms, jobname, options).run
      options = options.merge :mutated_isoforms => principal_isoforms, :database => database
      Structure.job(:annotated_variants, jobname, options)
    end
    task task_name => :tsv do 
      step(:annotated_variants).load
    end
    annotation_tasks << task_name

    task_name = [database,"neighbour_annotations"] * "_"
    task_name = task_name.to_sym

    dep do |jobname,options|
      principal_isoforms = GenomicTools::ExomeSequencing.job(:principal_isoforms, jobname, options).run
      options = options.merge :mutated_isoforms => principal_isoforms, :database => database
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
end
