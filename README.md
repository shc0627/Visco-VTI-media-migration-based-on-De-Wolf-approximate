####Application Name###
One-way Wave Migration Method for Visco-acoustic VTI Medium Based on De Wolf Approximation


#####Program Introduction######
The real earth medium is characterized by absorption attenuation and anisotropy, the overlook of which will lead to deviation in the amplitude and imaging position of the migration profile. There will be a large calculation error using the classical Born approximation method under the strong forward accumulative effect of a visco-acoustic VTI medium. Therefore, this paper first combines the De Wolf approximation with the visco-acoustic VTI medium and renormalizes the Born series by using the De Wolf approximation, thereby establishing the De Wolf approximate integral representation of the visco-acoustic VTI medium to improve the accuracy of the Born approximation. Then, the velocity model is divided into several thin slabs according to the theory of thin slab division. Medium parameters in the thin slabs are composed of background parameters (background velocity, background quality factor, and background anisotropy parameters) and disturbance parameters (velocity disturbance, quality factor, and anisotropy parameter disturbance). Afterwards, dual-domain screen approximation migration operators (frequency-wavenumber domain and frequency-space domain) of the visco-acoustic VTI medium are established. Finally, the prestack migration imaging method of visco-acoustic VTI medium based on De Wolf approximation is achieved.


### Program language ###
C Programming Language


#####Operating environment#####
Linux system, Ubuntu 18.04 or above


######Install OpenMPI######
Download URL: https://www.open-mpi.org/

(1)Upload the file to the server and decompress it

(2)Configure installation path, compile and install, customize installation path
./configure --prefix=/usr/local/openmpi
make
make install

(3)Set the environment variable to the path of your own installation
MPI_HOME=/usr/local/openmpi
export PATH=${MPI_HOME}/bin:$PATH
export LD_LIBRARY_PATH=${MPI_HOME}/lib:$LD_LIBRARY_PATH
export MANPATH=${MPI_HOME}/share/man:$MANPATH


###Program composition####
main: VTI_visco-Dewolf.c
subfunction: alloc.c, complex.c, pfafft.c
parameter file: vel.txt, epsilu.txt, deta.txt, factor.txt, record.txt
Parameter Table: pars.txt

     nz=220           Number of vertical grid points
     nx=600           Horizontal grid points
     dz=8             Vertical grid spacing
     dx=12            Horizontal grid spacing
     sx=300             Lateral coordinates of the earthquake source
     sz=0             Horizontal grid spacing
     nt=1701          Sampling length
     dt=0.001        Sampling interval
     fmax=20         Wavelet main frequency
     pxl=100          Expand the left model edge
     pxr=100          Expand the edge of the model on the right
     maxtrace=600    All earthquake record lengths
     ns=1             Total number of seismic sources
     ds=15             Cannon spacing
     np=1              Number of threads

output file: image_re_Dewolf.txt

####Program compilation and execution
Compile command
       gcc -O3 VTI_visco-Dewolf.c -o DVTI_visco-Dewolf.o -lm -g

Run command
       ./Dewolf.o


### Contact ###
Huachao Sun, shc0627@126.com
