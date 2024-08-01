# Correlation Feature Selection implementato in Linguaggio Assembly x86-32+SSE, x86-64+AVX e OpenMP

Il progetto ha come obiettivo quello di fornire una versione ottimizzata per il problema di "Correlation Feature Selection (CFS)".
L'ottimizzazione riguarda l'utilizzo di tecniche e meccanismi discussi nel file "relazione_progetto_architettute.pdf".

In generale, sono state presentate le seguenti implementazioni:
1. cfs32c.c - versione 32 bit del meccanismo di CFS
2. cfs64c.c - versione 64 bit del meccanismo di CFS
3. cfs32c_omp.c e cfs64c_omp.c - versione a 32 e 64 bit del meccanismo di CFS in cui si fa uso della libraria OMP.

L'ottimizzazione ha previsto la definizione di procedure in linguaggio assembly per l'esecuzione delle operazioni "dominanti" 
dell'algoritmo di CFS.

La repository contiene la cartella "dataset" in cui sono presenti dataset di test per CFS.
