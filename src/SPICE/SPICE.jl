# Instruments using SPICE save data into "experiment" directories, making figuring out where a scan is stored a bit easier
function SPICEdataPath(inst::AbstractString,experiment::Integer,no::Integer,base="/data/HFIR")
  # somehow have the default SPICE root stored in an environment variable?
  p=f=uppercase(inst)
  ("HB2C"==p || "WAND"==p) && (p="HB2C";w="WAND")
  joinpath(base,p,"exp$experiment","Datafiles",w*"_exp"*pad(experiment,4)*"_scan"*pad(no,4)*".dat")
end
function SPICEload(instrument::AbstractString,o...)
    instrument=uppercase(instrument)
    (ismatch(r"^HB2C",instrument)||ismatch(r"^WAND",instrument))&&(WANDload(o...))
end

include("SPICEload.jl") # a general SPICE file loader

include("WANDd.jl")
