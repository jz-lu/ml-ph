using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--deg", "-d"
            help = "twist angle in degrees"
            required = true
            arg_type = Float64
        "-N", "-n"
            help = "grid size (N if grid is NxN)"
            arg_type = Int
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    println(parsed_args["N"])
    println(parsed_args["deg"])
end

main()