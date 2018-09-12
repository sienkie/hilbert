# Hilbert Data Visualization

Galaxy tool visualization of NGS data on Hilbert Curve.

**What does it do?**
This tool maps genomic data stored in .bed or .bedgraph files on a Hilbert space-filling curve and creates an image file with visualization of the process. [Mapping algorithm](https://en.wikipedia.org/wiki/Hilbert_curve) provides a method to transfer between one-dimensional and two-dimensional space with preservation of locality. Taking that into account this tool allows the user to:
* inspect distribution of score throughout the whole genome,
* compare distribution of score between up to three datasets,
* easily inspect the quality of data.

The interface for this tool contains several options:

* **File with chromosome sizes** - Choose a .tsv file from history containing chromosomes sizes. Order of chromosomes in the file will be taken into account during visualization.
* **Input files** - Choose .bed or .bedgraph file from history to be visualized. You can specify a **color to plot dataset with** and the **threshold above which regions are colored**. When specified, regions with a score below the chosen threshold will be ignored during visualization. Use **Insert Input Files** option to add up to three datasets with separate colors and/or thresholds.
* **Create color map with chromosome territories** - If selected, the tool will create a separate .png file with chromosome territories map according to mapping algorithm.
* **Create text file with color legend** - If selected, the tool will create a separate .txt file with information about colors used to plot datasets with.
* **Invert colors** - If selected, the tool will generate printing friendly visualization on a white background instead of black. Colors used to plot datasets with also will be inverted. It is recommended to also use **create a text file with color legend** option to keep track of color changes.
* **Order of Hilbert Curve** - This value alongside reference chromosome sizes affect visualization precision. Using default value (11) results with creation of image with resolution 2048 x 2048px (2^11 x 2^11px). Every pixel corresponds to 2^11 x 2^11/(whole genome length) base pairs.


#### Installation

Currently package isn't available in the Galaxy Tool Shed. To install tools on individual instances of the Galaxy platform perform following steps:

* clone repository into your Galaxy tools/ directory:

~~~
git clone https://github.com/sienkie/hilbert.git
~~~

* copy lines below into your tool configuration file *tool_conf.xml* file located in the config/ directory of the Galaxy installation

~~~
<section id='ngsdatavis' name="NGS Data Visualization">
    <tool file="hilbert/hilbert.xml" />
</section>
~~~

For more information about configuring custom tools please refer to [this](https://galaxyproject.org/admin/tools/add-tool-tutorial/#4) tutorial.