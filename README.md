# ptx-multiband
This is a matlab toolbox for designing parallel transmit (pTx) multiband (MB) RF pulses.
You may run the demo script, `ptxMB.m`, to grab an idea of how this toolbox can be used to design band-specific pTx MB multispokes pulses by solving regularized minimization. 
You may run the demo script, `ptxMB_peakpowctrl.m`, to grab an idea of how this toolbox can be used to design band-specific pTx MB multispokes pulses by solving minimization with explicit peak RF power constraints (formulated based on composite multiband RF waveforms). 
Note that for the pulse designers included to run successfully, you will need to install the matlab software for disciplined convex programming (which can be found at http://cvxr.com/cvx/).

If you use the toolbox, please consider citing the following paper:
- Wu X, Schmitter S, Auerbach EJ, Moeller S, Ugurbil K, Van de Moortele PF. Simultaneous multislice multiband parallel radiofrequency excitation with independent slice-specific transmit B1 homogenization. Magn Reson Med 2013;70(3):630–638.

If you find the recipe for construction of peak power constraints helpful, please consider citing the following ISMRM abstract: 
- Wu, X., et al. Peak RF power constrained pulse design for multi-band parallel excitation, ISMRM 2013, p4253

### Copyright & License Notice
This software is copyrighted by Regents of the University of Minnesota and covered by US 10,684,337. Regents of the University of Minnesota will license the use of this software solely for educational and research purposes by non-profit institutions and US government agencies only. For other proposed uses, contact umotc@umn.edu. The software may not be sold or redistributed without prior approval. One may make copies of the software for their use provided that the copies, are not sold or distributed, are used under the same terms and conditions. As unestablished research software, this code is provided on an "as is'' basis without warranty of any kind, either expressed or implied. The downloading, or executing any part of this software constitutes an implicit agreement to these terms. These terms and conditions are subject to change at any time without prior notice.
