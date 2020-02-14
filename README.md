In both scripts define paths to signal and background files.
Both files must be normalized to effective luminocity=1.
For this moment scale factor 2.9 for background is used.
First step is to obtain boundaries for MVA. Then it's necessary
to insert them into script which optimize boundaries on MX.
To define minimal number of background events per category set minB.
This important for optimization of boundaries on MX. The definition
of initial values for MX is important. If the result of minimization
is shifted at more that 35 it's better to put initial values closer to
the estimated region and repeat procedure.

To optain boundaries for MVA: AnalCatNewMVA2A.C <br />
To optain boundaries for MX: AnalCatNewMX1.C <br />



