process INSTALL_DEPENDENCIES {
    /*
     * This process installs BLAST, and other dependencies,
     * so they are available for downstream processes.
     */
    
    tag "installing_dependencies"

    //conda 'bioconda::blast bioconda::mafft bioconda::trimal bioconda::raxml'
    // for container image issue -profile gcb or conda - TBD
    //conda "./metaphlan_env.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1':
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"


    output:
    path 'install_complete.txt'       , emit: dependenciesinstall
    //path 'dependencies_install_complete.txt', emit: dependenciesinstall


    script: 
    def args = task.ext.args ?: ''


    """
    echo "Install BLAST/ BLAST+ Suite for tools such as 'makeblastdb', 'blastn' etc."

    # Get architecture and OS
    ARCH=\$(uname -m)
    OS=\$(uname -s | tr '[:upper:]' '[:lower:]')

    # Normalise architecture naming to match filenames in BLAST directory
    if [ "\$ARCH" == "x86_64" ]; then
        ARCH="x64"
    fi

    # Fetch latest version
    BASE_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/"

    echo "Use OS, \${OS}, and architecture, \${ARCH}, to find appropriate file in \${BASE_URL}."
    
    # Silently list all files from BASE_URL, search and select most appropriate file for user's system using grep for basic regex matching; added extra escape '\' to interpret as bash
    LATEST_VERSION=\$(curl -s "\$BASE_URL" | grep -Eo "ncbi-blast-2\\.[0-9]+\\.[0-9]+\\+-\${ARCH}-\${OS}\\.tar\\.gz" | sort -V | tail -n 1)
    
    # Error message
    if [ -z "\$LATEST_VERSION" ]; then
        echo "Could not find the latest version of BLAST+."
        exit 1
    else
        echo "Selected file: \${LATEST_VERSION}"
    fi

    # Construct download URL, download and extract BLAST
    DOWNLOAD_URL="\${BASE_URL}\$LATEST_VERSION"
    echo "Selected file being downloading from: \$DOWNLOAD_URL"

    curl -L -O "\$DOWNLOAD_URL"
    tar -vxzf "\$LATEST_VERSION"

    # Use wildcard to match the extracted directory (assumes the prefix "ncbi-blast-" is consistent)
    EXTRACTED_DIR=\$(find . -maxdepth 1 -type d -name "ncbi-blast-*" | head -1)

    # Check if directory was extracted
    if [ -z "\$EXTRACTED_DIR" ]; then
        echo "Extracted directory not found."
        exit 1
    else
        # Rename the extracted directory to 'ncbi-blast'
        mv "\$EXTRACTED_DIR" ncbi-blast
    fi

    # Check existence of bin directory
    if [ ! -d "ncbi-blast/bin" ]; then
        echo "BLAST+ bin directory not found in renamed path."
        exit 1
    fi

    # Copy binaries to /usr/local/bin and clean up
    cp -r "ncbi-blast/bin/" "/usr/local/bin/"
    rm -rf "\$LATEST_VERSION"

    # Mark process completion
    touch install_complete.txt

    """
}

