import sys
def correct_spaces(ids_list)-> list:
    ids_list = [element.strip() for element in ids_list]
    ids_list[1] = " " + ids_list[1]
    return(ids_list)

def main():
    gtf_file = sys.argv[1]
    with open(gtf_file, "r") as f_in, open(gtf_file+".tmp", "w") as f_out:
        for line in f_in:
            ll = line.strip().split("\t")
            ids = ll[8]
            ids_list = ids.split(";")[0:2]
            if ids_list[0].startswith("transcript"):
                ids_list.reverse()
                ids_list = correct_spaces(ids_list)
                ll[8] = ";".join(ids_list)+";"
            new_line = "\t".join(ll)+"\n"
            f_out.write(new_line)

if __name__=="__main__":
    main()
