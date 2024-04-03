#!/usr/bin/env python

with open(species_file, "r") as f:
    for line in f:
        species, acc = line.rstrip().split("\t")
        if species == "dome":
            dome.append(acc)
        elif species == "wild":
            wild.append(acc)
        elif species == "ance":
            ance.append(acc)
        else:
            raise ValueError("Illegal species in species file. Can only be wild")

comp_species_dict = {"DOME": dome, "WILD": wild, "ANCE": ance}

chr_lens = dict()
with open(fai_file, "r") as f:
    for line in f:
        header, length = line.rstrip().split("\t")
        chr_lens[header] = int(length)

n_windows = chr_lens[chr_used] // window_size
comp_max_accs = np.zeros((len(dome), n_windows), dtype=np.int32)
current_window = 0
supp_matrix = []

with open(vcf_file, "r") as f:
    for line in f:
        line = line.rstrip()
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                acc_vec = [i.replace(".ont.s", "") for i in line.split("\t")[9:]]
                assert set(acc_vec) == set(dome + wild + ance)
        else:
            fields = line.split("\t")
            at_chr = fields[0]
                for j in tags:
                    if j.startswith("SUPP_VEC="):
                        supp_vec = [int(i) for i in j[9:]]
                if supp_vec is None:
                    raise ValueError("Missing 'SUPP_VEC' field")
                # Build the support matrix for this window or start a new matrix for a new window
                if widx == current_window:
                    supp_matrix.append(supp_vec)
                else:
                    sm = np.asarray(supp_matrix)
                    
                        for comp_acc in comp_species_dict[comp_species]:
                            comp_supp_idx = acc_vec.index(comp_acc)
                            this_comp_vec = sm[:, comp_supp_idx]
                            if np.count_nonzero(this_comp_vec) >= min_den and np.count_nonzero(this_vec) >= min_den:
                                num = np.count_nonzero(np.logical_and(this_vec, this_comp_vec))  # Intersection
                                den = np.count_nonzero(np.logical_or(this_vec, this_comp_vec))  # Union
                                t_distances[i].append(num / den)
                            else:
                                t_distances[i].append(-1)
                    t_distances_argmax = np.asarray([np.argmax(i) for i in t_distances])
                    comp_max_accs[:, current_window] = t_distances_argmax
                    t_distances = np.asarray([np.max(i) for i in t_distances])
                    distances[:, current_window] = t_distances
                    if widx == n_windows:
                        break
                    current_window = widx
                    supp_matrix = []
                    supp_matrix.append(supp_vec)

with open(simi_sp_tsv + "." + comp_species + ".tsv", "w") as f:
    f.write("Sample\t" + "\t".join( [str(i*window_size) for i in range(n_windows)]) + "\n")
    for i in range(len(dome)):
        f.write(dome[i] + "\t" + "\t".join( [comp_species_dict[comp_species][j] for j in list(comp_max_accs[i, :])] ) + "\n")

print("Sample\t" + "\t".join( [str(i*window_size) for i in range(n_windows)] ))
for i in range(len(dome)):
    print(dome[i] + "\t" + "\t".join( [str(j) for j in list(distances[i, :])] ).replace("-1.0", "NA"))

