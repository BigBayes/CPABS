
# Generate a smarter process manager, only spawns processes when there is an idle worker
function initialize_smart_spawn()

    jobs = Any[]
    running_procs = Dict{Int64, RemoteRef}()

    available_workers = workers()

    job_id_counter = 0

    function smart_spawn(args...)

        while length(available_workers) == 0
            sleep(100)
            for i in workers()
                if isready(running_procs[i])
                    push!(available_workers, i)
                end
            end
        end
        s = RemoteRef()
        worker = pop!(available_workers)
        @async put!(s, remotecall_fetch(worker, args...) )
        
        job_id_counter += 1
        push!(jobs, (worker, s))
        running_procs[worker] = s
        return job_id_counter, s
    end

    function get_job_counter()
        return job_id_counter
    end

    function get_jobs()
        return jobs
    end

    smart_spawn, get_job_counter, get_jobs
end

