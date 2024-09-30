process INSTALL_DEPENDENCIES {
    /*
     * This process installs BLAST, and other dependencies,
     * so they are available for downstream processes.
     */
    
    tag "installing_dependencies"

    output:
    path 'install_complete.txt'

    script: 

    """
    printf "Install BLAST/ BLAST+ Suite for tools such as 'makeblastdb', 'blastn' etc."
    curl -L -O "https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.15.0+-\$(uname)-\$(uname -m).tar.gz"
    tar -vxzf ncbi-blast-2.15.0+-\$(uname)-\$(uname -m).tar.gz
    sudo rm -rf ncbi-blast-2.15.0+-\$(uname)-\$(uname -m).tar.gz
    export PATH=\$PATH:\$PWD"/ncbi-blast-2.15.0+/bin"
    sudo cp "\$PWD/ncbi-blast-2.15.0+/bin/*" "/usr/local/bin/"

    printf "Installing other dependencies, mafft, trimal, raxml etc..."
    yes | conda install bioconda::mafft
    yes | conda install bioconda::trimal
    yes | conda install bioconda::raxml

    """
}
