"""
    overprint(str...; nblines=1)  

Print `str` over previous line.
"""
function overprint(str...; nblines=1)  
    for i in 1:nblines-1
        overprint(str...)
    end
    print("\u1b[1F")
    #Moves cursor to beginning of the line n (default 1) lines up   
    print(str...)   #prints the new line
   print("\u1b[0K\n") 
   # clears  part of the line.
   #If n is 0 (or missing), clear from cursor to the end of the line. 
   #If n is 1, clear from cursor to beginning of the line. 
   #If n is 2, clear entire line. 
   #Cursor position does not change. 
end

"""Clear `n` previous lines in terminal."""
function clearnlines(n)
    print("\u1b[", "$n", "F", "\u1b[J")
end

"""
    display_progress(i::Number, total::Number; n=10, name="", margin=20)

Display a nice progress bar based on `i / total`. Can be prefaced with `name`.
`margin` is the size of the left part area for the message.
`n` is the size of the bar in number of characters.
"""
function display_progress(i::Number, total::Number; n=10, name="", margin=20)
    if i == 1
        println()
    end
    print("\u1b[1F")
    printstyled(repeat(" ", margin - length(name))..., name, " ", color=:green)
    print("[", repeat("=", convert(Int64, round(n*(i/total)))), (i == total ? "=" : ">"), repeat(" ", n - convert(Int64, round(n*(i/total)))),  "] $i/$total")
    print("\u1b[0K\n") 
end