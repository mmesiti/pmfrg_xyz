module CPUInfo

export get_cpu_info

# level 1

"""
    get_cpu_frequency_info() -> CPUFreqInfo

Query CPU frequency information using `likwid-setFrequencies -p`.

Returns a Dict containing:
- Maximum frequency of HWThread 0
- Governor setting
- Turbo state
- Maximum uncore frequency for Socket 0
- name of the host

Assumes all HWThreads have the same configuration.

# Throws
- `ErrorException`: If `likwid-setFrequencies` command fails or output cannot be parsed
"""
function get_cpu_info()::Dict{String,Any}
    io = IOBuffer()
    run(pipeline(`likwid-setFrequencies -p`, stdout=io, stderr=io))
    output_str = String(take!(io))
    lines = split(output_str, '\n')

    hwp_on_line = find_line_with_pattern(lines,r"HWP")
    hwthread0_line = find_line_with_pattern(lines, r"^HWThread 0:")
    socket0_line = find_line_with_pattern(lines, r"^Socket 0:")

    governor = extract_governor(hwthread0_line)
    max_freq =extract_max_frequency(hwthread0_line)
    turbo = extract_turbo_state(hwthread0_line)


    uncore_max_freq = extract_uncore_max_frequency(socket0_line)

    hostname = strip(read(`hostname`,String))

    Dict("governor" => governor,
         "max_freq" => max_freq,
         "turbo" => turbo,
         "uncore_max_freq" => uncore_max_freq,
         "hostname" => hostname,
         "hwp on" => if hwp_on_line == "" "no" else "yes" end)
end

# level 2
function find_line_with_pattern(lines::Vector{SubString{String}}, pattern::Regex)::SubString{String}
    for line in lines
        if occursin(pattern, line)
            return line
        end
    end
    ""
end

function extract_max_frequency(hwthread_line::SubString{String})::Float64
    # Example: "HWThread 0: governor  performance min/cur/max 0.4/0.660316/2.6 GHz Turbo 1"
    # Extract the third number in min/cur/max pattern
    freq_match = match(r"min/cur/max\s+[\d.]+/([\d.]+)/([\d.]+)\s+GHz", hwthread_line)
    freq_match === nothing && error("Could not parse max frequency from: $hwthread_line")
    return parse(Float64, freq_match.captures[2])
end

function extract_governor(hwthread_line::SubString{String})::String
    # Example: "HWThread 0: governor  performance min/cur/max ..."
    # Extract governor name (word after "governor")
    gov_match = match(r"governor\s+(\w+)", hwthread_line)
    gov_match === nothing && error("Could not parse governor from: $hwthread_line")
    return gov_match.captures[1]
end

function extract_turbo_state(hwthread_line::SubString{String})::Bool
    # Example: "... GHz Turbo 1" or "... GHz Turbo 0"
    # Extract turbo state at end of line
    turbo_match = match(r"Turbo\s+([01])", hwthread_line)
    turbo_match === nothing && error("Could not parse turbo state from: $hwthread_line")
    return turbo_match.captures[1] == "1"
end

function extract_uncore_max_frequency(socket_line::SubString{String})::Float64
    # Example: "Socket 0: min/max 0.4/4.4 GHz"
    # Extract the second number in min/max pattern
    freq_match = match(r"min/max\s+([\d.]+)/([\d.]+)\s+GHz", socket_line)
    freq_match === nothing && error("Could not parse uncore max frequency from: $socket_line")
    return parse(Float64, freq_match.captures[2])
end

end # module CPUInfo
