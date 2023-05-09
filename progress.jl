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
function clearnlines(n)
    print("\u1b[", "$n", "F", "\u1b[J")
end

function display_progress(i, total; n=10)
    if i == 1
        println()
    end
    overprint("[", ["=" for j in 1:round(n*(i/total))]..., (i == total ? "=" : ">"), [" " for j in 1:n - round(n*(i/total))]...,  "] ", round(100*i/total, digits=2), " %")
end