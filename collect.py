import os,sys, shutil

if not os.path.exists("Input/Dataset"):
    print("Dataset not available")
    sys.exit(1)

datasetname = sys.argv[1]
datasetname = datasetname.strip()

print(datasetname)
print("Copying files...wait")

if datasetname != "Test30_CASP15" and datasetname != "ARES_benchmark2":
    print("Wrong dataset name, please enter Test30_CASP15 or ARES_benchmark2")
    sys.exit(1)

idlist = []

idpath = "Input/input.txt"

with open(idpath, "w") as wfile:

    if datasetname == "Test30_CASP15":

        datasetpath1 = "Input/Dataset/Test/Test30"

        for idname in os.listdir(datasetpath1):
            if len(idname) <= 6:
                idlist.append(idname)

        datasetpath2 = "Input/Dataset/Test/CASP15"

        for idname in os.listdir(datasetpath2):
            if len(idname) <= 6:
                idlist.append(idname)      

        counter = 0

        for idname in idlist:

            if counter <=29:
                decoypath = f"{datasetpath1}/{idname}"
            else:
                decoypath = f"{datasetpath2}/{idname}"

            for decoyname in os.listdir(decoypath):

                fulldecoypath = f"{decoypath}/{decoyname}"

                #print(fulldecoypath)

                dstpath = f"Input/RNA_pdbs/{idname}_{decoyname}"

                #print(dstpath)

                decoy_pfx = f"{idname}_{decoyname}".split(".")[0]

                #print(decoy_pfx)

                shutil.copy(fulldecoypath, dstpath)

                wfile.write(decoy_pfx + "\n")

                #break
            
            counter += 1


    else:

        datasetpath = "Input/Dataset/Ares_set/Ares_benchmark2/"

        for idname in os.listdir(datasetpath):
            if len(idname) == 4:
                idlist.append(idname)

        print(len(idlist))
        
        for idname in idlist:

            decoypath = f"{datasetpath}/{idname}"
            
            for decoyname in os.listdir(decoypath):

                fulldecoypath = f"{decoypath}/{decoyname}"

                #print(fulldecoypath)

                dstpath = f"Input/RNA_pdbs/{idname}_{decoyname}"

                #print(dstpath)

                decoy_pfx = f"{idname}_{decoyname}".split(".")[0]

                #print(decoy_pfx)

                shutil.copy(fulldecoypath, dstpath)

                wfile.write(decoy_pfx + "\n")

                #break

print("Done")