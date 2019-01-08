import sys
print(sys.argv)

if len(sys.argv) == 3:
    infn = sys.argv[1] 
    outfn = sys.argv[2] 
else:
    infn = "file.txt"
    outfn = "file.CDS.txt"

correct_header = True
c = w = 0

with open(outfn, "w") as f_out:
    for line in open(infn):
        if line.startswith(">"):
            header_fields = line.split("|")
            if len(header_fields)>=4:
                correct_header = True
                transcript_ID = header_fields[1]
                CDS_start = min([int(x) for x in header_fields[3].split(";")])
                CDS_end = max([int(x) for x in header_fields[2].split(";")])
                header = ">%s|%s|%s\n"%(transcript_ID, CDS_start, CDS_end)
                f_out.write(header)
                c += 1 
            else:
                correct_header = False
                w += 1
        else:
            if correct_header == True:
                f_out.write(line)
print(c, w)
