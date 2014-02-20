module GenomicTools
  def self.sample_info
    @sample_info ||= Persist.persist("Sample info", :marshal, :dir => persistdir) do
      samples = {}
      datadir = self.workdir
      datadir.glob("**/Samples/*").each do |file|
        relative_path = Misc.path_relative_to(datadir, file)
        center, experiment_type, file_type, name = relative_path.split("/")
        info = file.yaml
        info["CENTER_NAME"] = center
      samples[name] = IndiferentHash.setup info
      end
      samples
    end
  end

  def self.all_workflow_sample_files(workflow = ExomeSequencing)
    sample_files = {}
    workflow.tasks.each do |task_name,task|
      task_dir = workflow.workdir[task_name]
      job_names = workflow.jobs(task_name)
      job_files = Set.new

      job_names.each do |jn| 
        job = workflow.load_name(task_name, jn)
        job_files << job.path << job.info_file << job.files_dir 
      end.any?

      files = Dir.glob(File.join(task_dir, '*'))
      Log.debug("Found #{files.length} for #{ workflow } #{ task_name }")

      files.delete_if do |file|
        job_files.include? file
      end

      files.each do |file| 
        sample = File.basename(file)

        sample = (task.extension and not task.extension.empty?) ?
          sample.sub(/\.#{task.extension}$/,'') :
          sample

        sample_files[sample] ||= []
        sample_files[sample] << task_name 
      end
    end

    sample_files
  end

  def self.all_workflow_sample_jobs(workflow = ExomeSequencing)
    sample_jobs = {}
    workflow.tasks.each do |task_name,task|
      workflow.jobs(task_name).each do |sample|
        sample_jobs[sample] ||= []
        sample_jobs[sample] << task_name
      end
    end
    sample_jobs
  end
end
