Project developed as a Danuser Lab (Lyda Hill Department of Bioinformatics) research intern during 2025 UTSW STARS Summer Research Program. Created post acquisition registration routines to realign in-vivo selective plane illumination microscopy data of transgenic zerbrafish. Implemented a post acquisition  registration routine algorithm for realigning data using a temporal search and sort algorithm that I applied to create 2D and 3D timelapse visualizations for different data stacks acquired. 

**How to run**: Download code and place in the same directory as the selective plane illumination microscopy data. Note, this code is specifically tailored to the selective plane illumination microscopy data acquired and stored by the Danuser Lab and cannot be rerun unless in their environment and VPN. 

**Full Project Research Board**: https://drive.google.com/file/d/1NusKFMVcFF5HnGOu4XOgLQdVwwBy7jBy/view?usp=sharing

**Abstract**
Zebrafish can be used as a model for investigating the cardiovascular systems of human beings,
which is of interest for understanding the physiological origins of arrhythmias and other heart
diseases. Selective plane illumination microscopy (SPIM) facilitates fast imaging, low
photo-toxicity, low photobleaching, and good optical sectioning capability. However, since the
beating rate of zebrafish is about 2Hz, SPIM is still not fast enough to capture the zebrafish heart
in 3D before the heart expands and contracts, resulting in significantly misaligned raw images.
This study focuses on developing a pipeline to reconstruct images to visualize an aligned heartbeat
in 4D (xyzt). We use Pearson correlation, or Structural similarity index measure (SSIM), to realign
the images. We also apply the sliding window method to improve the reconstruction accuracy. We
hope our pipeline can help to reconstruct any zebrafish heart images acquired with SPIM, then we
can further perform cell segmentation and cell tracking, potentially also with machine learning
methods, to investigate the cardiovascular cycle. As a long-term goal, this work can help to study
cardiovascular irregularities or drug discovery for heart diseases.
