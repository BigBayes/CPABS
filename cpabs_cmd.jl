# run_phylo_cmd
# 
# Main entry point for running betasplitting phylogenetic reconstruction from the shell
#
#
 
using CPABS
using ArgParse

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin
        "input_file"
            help = "Input file containing read data (see README for format)" 
            required=true
        "--num-iterations", "-n"
            help = "Number of MCMC iterations to run (default 10000)"
            arg_type = Int
            default = 10000
        "--output", "-o"
            help = "Output file to store the MCMC chain"
        "--json_output", "-j"
            help = "Output directory to store the JSON summary"
        "--alpha", "-a"
            help = "alpha parameter, default = 0 (recommended)"
            arg_type = Float64
            default = 0.0
        "--multilocus-filename"
            help = "Input file for multilocus reads (see README for format)"
        "--wang-landau-partition", "-w"
            help = "Boundaries defining the WL partition (default = [4,5,6,7,Inf])"
            default = "[4,5,6,7,Inf]"
        "--init-K", "-K"
            help = "Number of clones used in the intial state (default = 3)"
            arg_type = Int
            default = 3


    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    inputfile = parsed_args["input_file"]
    num_iterations = parsed_args["num-iterations"]    
    outputfile = parsed_args["output"]
    json_dir = parsed_args["json_output"]
    alpha = parsed_args["alpha"]
    multilocus_filename = parsed_args["multilocus-filename"]
    wl_partition = parsed_args["wang-landau-partition"]
    init_K = parsed_args["init-K"]



    s = split(wl_partition[2:end-1], ",")
    wl_partition = [float(i)::Float64 for i in s] 

    num_iterations = num_iterations == nothing ? 10000 : num_iterations

    CPABS.run_phylo_experiment(inputfile, alpha, multilocus_filename,
                               wl_K_boundaries=wl_partition,
                               num_iterations=num_iterations,
                               outputfile=outputfile,
                               json_output_directory=json_dir,
                               init_K=init_K)
                          
end

main()
