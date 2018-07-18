status(level::AbstractString,o...)=status(Symbol(level),o...)
status(level::Symbol,message::AbstractString)=status(level,STDOUT,message)
function status(level::Symbol,io::IO,message::AbstractString)
    # colors known by `print_with_color` --
    # :normal, :default, :bold, :black, :blue, :cyan,
    # :green, :light_black, :light_blue, :light_cyan, :light_green, :light_magenta,
    # :light_red, :light_yellow, :magenta, :nothing, :red, :white, or :yellow
    known_levels=[:info,:update,:hint,:time,:debug]
    colors=[:default,:light_cyan,:light_magenta,:light_blue,:light_green,:light_red]
    i=findfirst(level.==known_levels)+1 # to ensure we can index into colors
    print_with_color(colors[i],io,Dates.format(now(),"yyyy-mm-dd HH:MM:SS "))
    print_with_color(colors[i],io,uppercase("$level:"),bold=true)
    print_with_color(colors[i],io," "*message*"\n")
end
export status
