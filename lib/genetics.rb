require 'rbbt/workflow'

module GenomicTools::Genetics
  extend Workflow

  extension :yaml
  task :Pedigree => :yaml do
    raise "No pedigree: #{ name }"
  end

  dep :Pedigree
  task :zygosis => :tsv do 
    pedigree = step(:Pedigree).load
    pedigree.extend IndiferentHash
    father =  pedigree[:father]
    mother =  pedigree[:mother]
    child =  pedigree[:child]

    father_mutations = Set.new GenomicTools::ExomeSequencing.job(:GenomicMutations, father).run
    mother_mutations = Set.new GenomicTools::ExomeSequencing.job(:GenomicMutations, mother).run
    child_mutations = Set.new GenomicTools::ExomeSequencing.job(:GenomicMutations, child).run
    
    tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Father", "Mother", "Novel", "Homozygous"])
    child_mutations.each do |mutation|
      from_father = father_mutations.include? mutation
      from_mother = mother_mutations.include? mutation
      tsv[mutation] = [from_father, from_mother,
        !(from_father or from_mother),
        (from_father and from_mother)].collect{|b| b.to_s }
    end

    tsv
  end

  dep :zygosis
  task :all => :string do
    "DONE"
  end
end
