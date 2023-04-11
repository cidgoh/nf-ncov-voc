# nf-ncov-voc Parameters

> Input/Output parameters

Define where the pipeline should find input data and save output data.<br>

<ul>
<li><code> --startdate </code> Starting date to extractdataset
(yyyy-mm-dd). </li>
<li><code> --enddate </code> Starting date to extractdataset
(yyyy-mm-dd). </li>

</ul>

> Quality Control parameters (`BBMap`)

Define where the pipeline should find input data and save output data.<br>

<ul>
<li><code> --maxns </code> Reads shorter than this after trimming will
be discarded. </li>
<li><code> --minlength </code> Reads shorter than this after trimming
will be discarded. </li>
</ul>

> Mapping parameters (`MiniMap2`/`BWA`)

Define where the pipeline should find input data and save output data.<br>

<ul>
<li><code> --ref </code> Instead of indexing the reference file in the
again, the prefix of previously-created reference index files can be
used. </li>
<li><code> --keep_min_map_quality </code> Minimum mapping quality of
covid reads
to keep </li>
<li><code> --remove_min_map_quality </code> Minimum mapping quality
of the human reads to remove </li>

</ul>

> Variant Calling parameters (`Freebayes`/`iVar`)

Define where the pipeline should find input data and save output data.<br>

<ul>
<li><code> --ploidy </code> Sets the ploidy for the analysis </li>
<li><code> --var_MinFreqThreshold </code>  Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position. </li>
<li><code> --var_MinDepth </code> Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position. </li>
<li><code> --mpileupDepth </code> Mpileup depth for iVar (although undocumented in mpileup, setting to zero removes limit). </li>
<li><code> --var_FreqThreshold </code> Frequency threshold for consensus variant. </li>
<li><code> --var_MinVariantQuality </code> Minimum mapQ to call variant. </li>
<li><code> --lower_ambiguityFrequency </code> Variants with frequency less that this will be discarded. </li>
<li><code> --upper_ambiguityFrequency </code> Substitution variants with frequency less than this will be encoded with IUPAC ambiguity codes </li>
</ul>
