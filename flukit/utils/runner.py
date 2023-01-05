import shutil

def find_executable(names, default = None):
    """
    Return the path to the first executable found in PATH from the given list
    of names.
    Raises a (hopefully helpful) error if no executable is found.  Provide a
    value for the "default" parameter to instead return a value.
    """
    exe = next(filter(shutil.which, names), default)

    if exe is None:
        print("Unable to find any of %s in PATH=%s" % (names, os.environ["PATH"]))
        print("\nHint: You can install the missing program using conda or homebrew or apt-get.\n")
        raise Exception

    return exe

def build_nextclade():
    '''
    build nextclade arguments
    '''

def build_iqtree(
	aln_file, 
	out_file, 
	substitution_model="GTR", 
	clean_up=True, 
	nthreads=1, 
	tree_builder_args=None):
    '''
    build tree using IQ-Tree with parameters "-fast"
    arguments:
        aln_file    file name of input aligment
        out_file    file name to write tree to
    '''

    # find executable in path
    iqtree = find_executable([
        "iqtree2",
        "iqtree"
    ])
   
    # construct the command syntax
    if substitution_model.lower() != "auto":
        call = [iqtree, "-ntmax", str(nthreads), "-s", shquote(tmp_aln_file),
                "-m", substitution_model, tree_builder_args, ">", log_file]
    else:
        call = [iqtree, "-ntmax", str(nthreads), "-s", shquote(tmp_aln_file), tree_builder_args, ">", shquote(log_file)]

    cmd = " ".join(call)

    # print a heplful message
    print("Building a tree via:\n\t" + cmd +
          "\n\tNguyen et al: IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies."
          "\n\tMol. Biol. Evol., 32:268-274. https://doi.org/10.1093/molbev/msu300\n")
    if substitution_model.lower() == "auto":
        print(f"Conducting a model test... see '{shquote(log_file)}' for the result. You can specify this with --substitution-model in future runs.")

    # run
    try:
        run_shell_command(cmd, raise_errors = True)
        T = Phylo.read(tmp_aln_file+".treefile", 'newick')
        shutil.copyfile(tmp_aln_file+".treefile", out_file)
        for n in T.find_clades(terminal=True):
            tmp_name = n.name
            for v,c in reverse_escape_dict.items():
                tmp_name = tmp_name.replace(v,c)
            n.name = tmp_name
        #this allows the user to check intermediate output, as tree.nwk will be
        if clean_up:
            if os.path.isfile(tmp_aln_file):
                os.remove(tmp_aln_file)

            for ext in [".bionj",".ckp.gz",".iqtree",".mldist",".model.gz",".treefile",".uniqueseq.phy",".model"]:
                if os.path.isfile(tmp_aln_file + ext):
                    os.remove(tmp_aln_file + ext)
    except:
        print("ERROR: TREE BUILDING FAILED")
        if os.path.isfile(log_file):
            print("Please see the log file for more details: {}".format(log_file))
        T=None
    return T
