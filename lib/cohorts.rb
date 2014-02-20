
Workflow.require_workflow "Enrichment"
module GenomicTools::Cohorts
  extend Workflow

  task :Cohort => :array do
    raise "Cohort needs to be defined as a list of samples"
  end

  dep :Cohort
  task :GeneCounts => :tsv do
    cohort = step(:Cohort).load
    all_genes = []
    cohort.each do |sample|
      all_genes.concat GenomicTools.sample_job(GenomicTools::ExomeSequencing, sample, "AffectedGenes").load
    end
    TSV.setup(Misc.counts(all_genes), :key_field => "Ensembl Gene ID", :fields => ["Count"], :cast => :to_i, :type => :single)
  end

  dep :GeneCounts
  input :recurrence_threshold, :integer, "Minimun sample count for recurrence", 2
  task :RecurrentGenes => :array do |recurrence_threshold|
    counts = step(:GeneCounts).load
    counts.select("Count"){|c| c.to_i >= recurrence_threshold }.keys
  end

  dep :RecurrentGenes
  input :cutoff, :float, "Enrichment cutoff (FDR q-value)", 0.1
  task :PathwayEnrichment => :string do |cutoff|
    genes = step(:RecurrentGenes).load
    databases = %w(kegg go_bp)
    databases.each do |database|
      Open.write(file(database), Enrichment.job(:enrichment, name, :database => database, :cutoff => cutoff, :list => genes).run.to_s)
    end
    "Enrichment done"
  end

  dep :RecurrentGenes
  task :all => :string do
    "DONE"
  end

end
