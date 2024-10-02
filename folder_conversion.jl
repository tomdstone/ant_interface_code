using ArgParse

s = ArgParseSettings(
    description = "This program processes all `.cnt` files in a folder structure, converting them to `.set`"
    , version = "1.1"
)

@add_arg_table! s begin
    "dir"
        arg_type = String
        help = "Root directory to process all folders in"
        default = "."
    "-v", "--version"
        action = :show_version
        help = "show version information and exit"
end

dir = abspath(parse_args(s)["dir"])

for (r, _, files) in walkdir(dir)
    for file in files
        if splitext(file)[2] == ".cnt" && !(splitext(file)[1] * ".set" in files)
            infile  = joinpath(r, file)
            outfile = joinpath(r, splitext(file)[1]*".set")
            command = `matlab -nodisplay -nosplash -nodesktop -r "convert_cli $infile $outfile; exit;"`

            run(command)
        end
    end
end