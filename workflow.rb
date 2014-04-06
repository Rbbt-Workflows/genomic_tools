require 'rbbt/workflow'

Workflow.require_workflow "Genomics"
Workflow.require_workflow "Sequence"
Workflow.require_workflow "Structure"
Workflow.require_workflow "Appris"

module GenomicTools
  extend Workflow
end

require 'utils'
require 'exome_sequencing'
require 'cohorts'
require 'genetics'

module GenomicTools
  class << self
    attr_accessor :nested_workflows
  end
  self.nested_workflows = []

  def self.nest_workflow(workflow)
    nested_workflows << workflow

    class << workflow
      def workdir
        GenomicTools.workdir[self.to_s.split(":").last]
      end

      def workdir=(workdir)
        GenomicTools.workdir = workdir
      end
    end
  end

  self.nest_workflow ExomeSequencing
  self.nest_workflow Cohorts
  self.nest_workflow Genetics
end
