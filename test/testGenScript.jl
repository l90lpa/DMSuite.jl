using MATLAB

# Location of the Matlab implementation of the dmsuite library. Current the folder is expected to be in the same 
# directory as this script
MATLAB_DMSUITE_DIR = (@__DIR__)*"/dmsuite"

FFT_PACKAGE = "FFTW"



##### Helper Functions

function createTest(file, session, tolerance, outCount, funcname, args)
    eval_string(s,"clear")
    
    argCount = length(args)
    outNames = "out".*string.(1:outCount)
    argNames = "arg".*string.(1:argCount)

    outVars = computeMatlabResults(s,outNames,funcname,args,argNames)
    
    ### Write the test

    # serialise matlab results as test data
    testVarPrefix = "test_"
    for i in 1:outCount
        name = testVarPrefix * outNames[i]
        var = outVars[i]
        if length(size(var)) == 3
            for j in 1:size(var,3)
                writeVar(file, name * "_" * string(j), var[:,:,j])
                write(file,"\n")
            end
        else
            writeVar(file, name, var)
            write(file,"\n")
        end
    end
    write(file,"\n")
    
    # add argument variables
    for i in 1:argCount
        writeVar(file, argNames[i], args[i]);  write(file,"\n")
    end
    write(file,"\n")
    
    # add Julia function evaluation
    write(file,join(outNames,",")*"="*funcname*"("*join(argNames,",")*")\n")
    write(file,"\n")

    # split any output variable that is a three dimensional array (into multiple matrices)
    split = false
    for i in 1:outCount
        var = outVars[i]
        name = outNames[i]
        if length(size(var)) == 3
            for j in 1:size(var,3)
                write(file,name*"_"*string(j)*"="*name*"[:,:,"*string(j)*"]\n")
            end
            split = true
        end
    end
    if split; write(file,"\n") end
    
    # add tests for testing Julia results against matlab results
    tol = string(tolerance)
    for i in 1:outCount
        var = outVars[i]
        name = outNames[i]
        if length(size(var)) == 3
            for j in 1:size(var,3)
                suffix = "_"*string(j)
                writeCompTest(file,testVarPrefix*name*suffix,name*suffix,tol)
                write(file,"\n")
            end
        else
            writeCompTest(file,testVarPrefix*name,name,tol)
            write(file,"\n")
        end
    end
    write(file,"\n")
end

function computeMatlabResults(s,outNames,funcname,args,argNames)
    for i in 1:length(args)
        arg = args[i]
        if typeof(arg)==Int64
            arg = Float64(arg)
        end 
        put_variable(s, Symbol(argNames[i]), arg)
    end

    eval_string(s,"["*join(outNames,",")*"]="*funcname*"("*join(argNames,",")*");")

    # Extract the matlab results

    outVars = []
    for i in 1:length(outNames)
        r = jarray(get_mvariable(s, Symbol(outNames[i])))
        push!(outVars, r)
    end
    outVars
end

function writeCompTest(file,varname1,varname2,tol)
    write(file,"@test maximum(abs.(" * varname1 * " - " * varname2 * ")) <= " * string(tol))
end

function writeVar(file,name,value)
    write(file, name*"="*string(value))
end



##### Define the tests

tol = 1e-15;
s = MSession()
eval_string(s,"addpath('"*MATLAB_DMSUITE_DIR*"')")

file = open("runtests.jl","w")
write(file, "##### This file has been generated by running the file testGenScript.jl\n")
write(file, "using Test\n")
write(file, "using "*FFT_PACKAGE*"\n")
write(file, "using DMSuite\n\n\n")
write(file, "@testset \"DMSuite Tests\" begin\n\n")

createTest(file, s, tol, 2, "chebdif", (8,2))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 1, "chebdifft", ([1.0 2.0 3.0 4.0],3))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 1, "chebint", ([1.0 2.0 3.0 4.0 5.0 4.0],[-1 -0.5 0 0.5 1]))
write(file, "\n\n############\n\n")

N = 8
g = [1.0 0 1.0; 1.0 0 1.0]
createTest(file, s, tol, 5, "cheb2bc", (N,g))
write(file, "\n\n############\n\n")
g = [0 1.0 1.0; 1.0 0 1.0]
createTest(file, s, tol, 5, "cheb2bc", (N,g))
write(file, "\n\n############\n\n")
g = [1.0 0 1.0; 0 1.0 1.0]
createTest(file, s, tol, 5, "cheb2bc", (N,g))
write(file, "\n\n############\n\n")
g = [1.0 1.0 1.0; 1.0 1.0 1.0]
createTest(file, s, tol, 5, "cheb2bc", (N,g))
write(file, "\n\n############\n\n")

createTest(file, s, tol, 2, "cheb4c", (9))
write(file, "\n\n############\n\n")

createTest(file, s, tol, 2, "fourdif", (8,0))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 2, "fourdif", (8,1))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 2, "fourdif", (9,1))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 2, "fourdif", (8,2))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 2, "fourdif", (9,2))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 2, "fourdif", (8,3))
write(file, "\n\n############\n\n")
createTest(file, s, 5e-15, 2, "fourdif", (9,3))
write(file, "\n\n############\n\n")

createTest(file, s, tol, 1, "fourdifft", ([1.0 2.0 3.0 4.0 5.0 6.0],3))
write(file, "\n\n############\n\n")

createTest(file, s, tol, 1, "fourint", ([1.0 2.0 3.0 4.0 5.0 6.0],[0 0.25 0.5 0.75 1.0]))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 1, "fourint", ([1.0 2.0 3.0 4.0 5.0],[0 0.25 0.5 0.75 1.0]))
write(file, "\n\n############\n\n")

createTest(file, s, tol, 2, "sincdif", (9,2,0.5))
write(file, "\n\n############\n\n")
createTest(file, s, 5e-15, 1, "sincdifft", ([1.0 2.0 3.0 4.0 5.0],2,0.5))
write(file, "\n\n############\n\n")

createTest(file, s, tol, 1, "herroots", (8))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 1, "lagroots", (8))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 1, "legroots", (8))
write(file, "\n\n############\n\n")

createTest(file, s, tol, 2, "herdif", (8,2,0.7))
write(file, "\n\n############\n\n")
createTest(file, s, tol, 2, "lagdif", (8,2,0.7))
write(file, "\n\n############\n\n")


write(file, "end")
close(file)

close(s)