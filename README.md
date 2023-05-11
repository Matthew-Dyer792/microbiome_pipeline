# Microbiome Pipeline
We present a novel microbiome analysis pipeline that quantifies microbial species counts from raw sequencing reads not belonging to the target organism's genome. This pipeline leverages the R package MetaScope to further reduce the error in the final output counts. MetaScope is a tool that uses a probabilistic model to assign reads to microbial taxa based on their similarity to reference genomes. Our pipeline integrates MetaScope with other methods for quality control, data preprocessing, and network analysis, and provides a comprehensive framework for microbiome studies. We demonstrate the performance and applicability of our pipeline on several datasets from different environments and organisms.

## Nextflow
Nextflow is a workflow management system that allows us to write scalable and reproducible pipelines for analysis. It simplifies the execution of complex tasks across different computing platforms and environments. Singularity is a container manager that enables us to run applications in isolated and portable environments. It helps us to avoid dependency conflicts and ensure reproducibility of our results. Anaconda is a software distribution that provides us with a large collection of tools and libraries for data analysis, including microbiome-specific packages. It also allows us to create and manage virtual environments for different projects. Together, these three software's complement each other by providing us with a flexible, robust and efficient framework for microbiome analysis.

# Installation
## Anaconda
We suggest using anaconda as it is a distribution of packages built for data science, and comes with conda, a package, and environment manager that you can use to create environments for isolating your projects that use different versions of Python and/or different version of packages.

1. Download the Anaconda installer for Linux from the [official website](https://docs.anaconda.com/anaconda/install/index.html).

2. Open a terminal window and navigate to the directory where you downloaded the installer.

3. Run the following command to add executable permission to the installer:
    ```bash
    chmod +x Anaconda3-YOUR_VERSION-Linux-x86_64.sh
    ```

4. Run the installer by running the following command:
    ```bash
    ./Anaconda3-YOUR_VERSION-Linux-x86_64.sh
    ```

5. Follow the prompts on the installer screens.

6. Once installation is complete, you can start using Anaconda by opening a new terminal window.

## Nextflow
Before installing Nextflow, you will need to make sure that **Java version 11** or greater is installed on your machine. You can check your Java version by running the following command in your terminal window:
```bash
java -version
```

1. Download the Nextflow executable by copying and pasting one of the following commands in your terminal window:
    ```bash
    wget -qO- https://get.nextflow.io | bash
    ```
    or
    ```bash
    curl get.nextflow.io | bash
    ```

2. Make the binary executable on your system by running:
    ```bash
    chmod +x nextflow
    ```

3. Optionally, move the nextflow file to a directory accessible by your $PATH:
    ```bash
    sudo mv nextflow /usr/local/bin/
    ```

4. If you want to install Nextflow through Anaconda, you can do so by running:
    ```bash
    conda create -n nextflow -c bioconda nextflow
    ```

## Singularity
Singularity is a free and open-source computer program that performs operating-system-level virtualization also known as containerization. One of the main uses of Singularity is to bring containers and reproducibility to scientific computing and the high-performance computing (HPC) world1. Singularity is a container framework designed to run scientific applications on HPC-backed resources.

Singularity is chosen by most HPCs as their primary container software because it allows users to pack an application/workflow/pipeline and all of its dependencies into a single image (file), which can be easily transported between different HPC systems3. Furthermore, Singularity assumes that the user does not have root privileges on the host OS, which makes it more secure than other containerization technologies.

1. Install dependencies:
    ```bash
    sudo apt-get update
    sudo apt-get install -y build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools libseccomp-dev wget pkg-config git cryptsetup-bin
    ```

2. Install Go:
    - Download the Go binary for Linux by going to the site [go.dev](https://go.dev/) and then clicking on Download Go.

    - Extract the Golang binaries tarball using the tar command to a directory of your choice:
        ```bash
        tar -C /usr/local -xzf go$VERSION.$OS-$ARCH.tar.gz
        ```
    - Extract the Golang binaries tarball using the tar command to a directory of your choice:
        ```bash
        export PATH=$PATH:/usr/local/go/bin
        ```
    - Extract the Golang binaries tarball using the tar command to a directory of your choice:
        ```bash
        go version
        ```

3. Download Singularity from a release:
    ```bash
    VERSION=3.8.0
    wget https://github.com/hpcng/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
    tar -xzf singularity-${VERSION}.tar.gz
    cd singularity-${VERSION}
    ```

4. Compile Singularity:
    ```bash
    ./configure --prefix=/usr/local
    make
    sudo make install
    ```

5. Verify that Singularity is installed correctly:
    ```bash
    singularity --version
    ```

Alternatively if you are using CAIR you can enter the following commands to activate singularity:
```bash
module load application/go/1.14.2
module load application/singularity/3.5.3 
```
# Running the Pipeline

```bash
conda activate nextflow
nextflow run /research/project/shared/benoukraf_lab/nextflow/microbiome_pipeline/main.nf -profile chia,ont,conda --input fastq/00811_8.fastq --target_index "/research/project/shared/benoukraf_lab/pathoscope/metascope/index/refseq_all/ont/*-ont.mmi" --filter_index /research/project/shared/benoukraf_lab/pathoscope/metascope/index/human/ont/hg38-ont.mmi

/research/project/shared/benoukraf_lab/nextflow/microbiome_pipeline/main.nf
-profile chia,ont,conda
--input fastq/00811_8.fastq
--target_index "/research/project/shared/benoukraf_lab/pathoscope/metascope/index/refseq_all/ont/*-ont.mmi"
--filter_index /research/project/shared/benoukraf_lab/pathoscope/metascope/index/human/ont/hg38-ont.mmi
```