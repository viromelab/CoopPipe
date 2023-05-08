# CoopPipe
Cooperative Pipeline for human viral genome reconstruction.

### REPLICATION ###

To download CoopPipe in a Linux system, please run:
<pre>
git clone https://github.com/mirakaya/CoopPipe.git
cd CoopPipe/src/
chmod +x *.sh
</pre>

If Miniconda is not present you can download and install it with:
<pre>
./CoopPipe --miniconda
</pre>

Install all of the tools necessary to execute the pipeline with:
<pre>
./CoopPipe -i -t 4
</pre>

To execute the whole pipeline, type:
<pre>
./CoopPipe.sh -r1 reads_forward.fq -r2 reads_reverse.fq
     --output out_analysis --threads 4 --memory 28        
     --tools HVRS/src --classifier falcon --all
</pre>

### CITATION ###

On using this software/method please cite:

* pending

### ISSUES ###

For any issue let us know at [issues link](https://github.com/mirakaya/CoopPipe/issues).

### LICENSE ###

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>
