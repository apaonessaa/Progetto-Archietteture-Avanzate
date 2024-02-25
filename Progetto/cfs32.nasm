; ---------------------------------------------------------
; Regressione con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
	maskres 	dd 		0x0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF
	maskfloat 	dd		0x3f800000, 0x0, 0x0, 0x0	
	maskabs		dd 		0x7FFFFFFF, 0x0, 0x0, 0x0			

section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	sc			resd	1
	mem_rff 	resd	1

section .text			; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in eax
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global prova

input		equ		8

msg	db	'sc:',32,0
nl	db	10,0

prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------
		; elaborazione
		; esempio: stampa input->sc
		; mov eax, [ebp+input]	; indirizzo della struttura contenente i parametri
        ; [eax] input->ds; 			// dataset
		; [eax + 4] input->labels; 	// etichette
		; [eax + 8] input->out;	// vettore contenente risultato dim=(k+1)
		; [eax + 12] input->sc;		// score dell'insieme di features risultato
		; [eax + 16] input->k; 		// numero di features da estrarre
		; [eax + 20] input->N;		// numero di righe del dataset
		; [eax + 24] input->d;		// numero di colonne/feature del dataset
		; [eax + 28] input->display;
		; [eax + 32] input->silent;
		; MOVSS xmm0, [eax+12]
		; MOVSS [sc], xmm0 
		; prints msg            
		; printss sc     
		; prints nl
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi			; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp			; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante

global mediaP4

mediaP4:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; media(fx, n%p, mu*, input*, strcf*)
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; Elaborazione
		; ------------------------------------------------------------
		mov eax, [ebp + 20]			; input*
		mov ebx, [eax + 20]			; input->N
		mov edx, [ebp + 12]			; N%p 
		neg edx						; -N%p 
		add edx, ebx				; N - N%p
		mov ecx, [ebp + 8]			; fx
		imul ebx, ecx				; N * fx
		mov ecx, [eax + 4]			; input->labels*
		mov eax, [eax]				; input->ds*
		xor esi, esi				; i=0
		xorps xmm5, xmm5			; n1
		xorps xmm6, xmm6			; m1
		xorps xmm7, xmm7			; m0
		xorps xmm3, xmm3			; 0 0 0 0

for_media:
		cmp esi, edx				; i < N-N%p
		jnl	res_media
		movups xmm0, [eax + ebx*4]	; [ ds[fx*N+i] | ds[fx*N+i+1] | ds[fx*N+i+2] | ds[fx*N+i+3] ]
		movups xmm1, xmm0			; [ ds[fx*N+i] | ds[fx*N+i+1] | ds[fx*N+i+2] | ds[fx*N+i+3] ]
		movups xmm2, [ecx + esi*4]	; [ labels[i]  | labels[i+1]  | labels[i+2]  | labels[i+3]  ]

		addps xmm5, xmm2			; counter n1 
		cmpneqps xmm2, xmm3			; labels_mask based on c=1 

		andps xmm0, xmm2			; ds and labels_mask	  (c=1)
		andnps xmm2, xmm1			; ds and not labels_mask  (c=0)

		addps xmm6, xmm0			; m1 m1 m1 m1
		addps xmm7, xmm2			; m0 m0 m0 m0

		add esi, 4					; i+=4
		add ebx, 4					; (fx*N+i)+4
		jmp for_media

res_media:
		mov edx, [ebp + 12]			; N%p
		cmp edx, 0							
		je end_media
		neg edx 					; -N%p
		add edx, 4					; p-N%p

		sub ebx, edx 				; fx*N+i - (p-N%p)
		sub esi, edx				; i - (p-N%p)

		movups xmm0, [eax + ebx*4]	; [ ds[fx*N+(N-4)] | ds[fx*N+(N-3)] | ds[fx*N+(N-2)] | ds[fx*N+(N-1)] ]
		movups xmm1, xmm0			; [ ds[fx*N+(N-4)] | ds[fx*N+(N-3)] | ds[fx*N+(N-2)] | ds[fx*N+(N-1)] ]
		movups xmm2, [ecx + esi*4]	; [ labels[N-4]  | labels[N-3]  | labels[N-2]  | labels[N-1]  ]

		movups xmm3, [maskres]		; residue_mask: 1 1 1 0
		
		cmp edx, 2					; p-N%p == 2
		jne res1_media
		insertps xmm3, xmm3, 0x3	; residue_mask: 1 1 0 0
		jmp calc_res_media

res1_media:	
		cmp edx, 3					; p-N%p == 3
		jne calc_res_media
		insertps xmm3, xmm3, 0x7	; residue_mask: 1 0 0 0

calc_res_media:
		andps xmm0, xmm3			; residue_mask on [ ds[fx*N+(N-4)] | ds[fx*N+(N-3)] | ds[fx*N+(N-2)] | ds[fx*N+(N-1)] ]
		andps xmm1, xmm3			; residue_mask on [ ds[fx*N+(N-4)] | ds[fx*N+(N-3)] | ds[fx*N+(N-2)] | ds[fx*N+(N-1)] ]
		andps xmm2, xmm3			; residue_mask on [ labels[N-4]  | labels[N-3]  | labels[N-2]  | labels[N-1]  ]
		addps xmm5, xmm2			; counter n1

		xorps xmm4, xmm4			; 0 0 0 0
		cmpneqps xmm2, xmm4			; labels_mask based on c=1 
		andps xmm2, xmm3			; residue_mask for labels_mask

		andps xmm0, xmm2			; ds and labels_mask	  (c=1)
		andnps xmm2, xmm1			; ds and not labels_mask  (c=0)
		
		addps xmm6, xmm0			; m1 m1 m1 m1
		addps xmm7, xmm2			; m0 m0 m0 m0

end_media:	
		xorps xmm0, xmm0
		xorps xmm1, xmm1
		mov eax, [ebp + 20]			; input*

		haddps xmm5, xmm5			; x x n1 n1
		haddps xmm5, xmm5			; x x x n1

		haddps xmm6, xmm6			; x x m1 m1
		haddps xmm6, xmm6			; x x x m1

		haddps xmm7, xmm7			; x x m0 m0
		haddps xmm7, xmm7			; x x x m0

		movss xmm0, xmm6			; copy m1
		divss xmm0, xmm5			; x x x m1/n1=mu1

		cvtsi2ss xmm1, [eax + 20]	; input->N Convert_Integer_To_Single_Precision_Floating_Point
		movss xmm4, xmm1			; copy N

		addss xmm6, xmm7			; x x x m1+m0
		divss xmm6, xmm1			; x x x (m1+m0)/N=mu

		subss xmm1, xmm5 			; x x x N-n1(=n0)
		divss xmm7, xmm1			; x x x m0/n0=mu0

		mov eax, [eax + 24]			; input->d
		mov ecx, [ebp + 8]			; fx
		mov ebx, [ebp + 16]			; mu*

		movss [ebx + ecx*4], xmm7	; mu[fx]=mu0
		add ecx, eax				; fx+d
		movss [ebx + ecx*4], xmm0	; mu[fx+d]=mu1
		add ecx, eax				; fx+d+d	
		movss [ebx + ecx*4], xmm6	; mu[fx+2*d]=mu

		subss xmm7, xmm0			; mu0-mu1
		movups xmm3, [maskabs]		; load mask abs
		andps xmm7, xmm3			; |mu0-mu1|
		divss xmm7, xmm4			; |mu0-mu1|/N

		sqrtss xmm1, xmm1			; sqrt(n0)
		sqrtss xmm5, xmm5			; sqrt(n1)

		subss xmm4, [maskfloat]		; x x x N-1
		sqrtss xmm4, xmm4			; sqrt(N-1)

		mulss xmm7, xmm1			; |mu0-mu1|/N * sqrt(n0)
		mulss xmm7, xmm5			; |mu0-mu1|/N * sqrt(n0) * sqrt(n1)
		mulss xmm7, xmm4			; |mu0-mu1|/N * sqrt(n0) * sqrt(n1) * sqrt(N-1)

		mov ecx, [ebp + 8]			; fx
		mov ebx, [ebp + 24]			; strcf*
		movss [ebx + ecx*4], xmm7	; strcf[fx] = |mu0-mu1|/N * sqrt(n0) * sqrt(n1) * sqrt(N-1)
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante

global rcfP4

rcfP4:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; rcfP4(d%2, N%4, mu, input, strcf)
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; Elaborazione
		; ------------------------------------------------------------
		mov ecx, [ebp+20]					; input
		mov eax, [ecx+20]					; N
		mov ecx, [ecx+24]					; d
		mov esi, [ebp+8]					; d%p
		neg esi 							; -d%p
		add esi, ecx						; d-d%p
		xorps xmm7, xmm7					; rcfMax
		push 0								; fmax
		xor edi, edi
for_rcf: 
		cmp edi, esi						; fi < d-d%p 
		jge resD_rcf
		mov edx, [ebp+16]					; mu*
		mov ecx, [ebp+20]					; input
		mov ecx, [ecx+24]					; d
		shl ecx, 1							; 2*d
		add ecx, edi						; 2*d+fi
		movss xmm0, [edx+ecx*4]				; x x x muu (fi)
		shufps xmm0, xmm0, 0				; muu[] (fi)
		movss xmm4, [edx+ecx*4+4]			; x x x muu (fi+1)
		shufps xmm4, xmm4, 0				; muu[] (fi+1)

		xorps xmm1, xmm1					; ximuu (fi)
		xorps xmm5, xmm5					; ximuu (fi+1)
		xor eax, eax						; contatore for innestato
		mov ecx, [ebp+20]					; input
		mov edx, [ecx+20]					; N
		mov ebx, [ebp+12]					; N%p
		neg ebx								; -N%p
		add ebx, edx						; N-N%p
		mov ecx, [ecx]						; ds
		mov esi, edi						; copy fi
		inc esi								; fi+1
		imul esi, edx						; (fi+1)*N
		imul edx, edi						; fi*N

forN_rcf: 
		cmp eax, ebx						; i < N-N%p
		jge resN_rcf

		movups xmm2, [ecx+edx*4]			; ds[fi*N+i...i+4]
		subps xmm2, xmm0					; x-muu
		mulps xmm2, xmm2 					; (x-muu)*(x-muu)
		addps xmm1, xmm2					; ximuu+=(x-muu)*(x-muu) (fi)

		movups xmm2, [ecx+esi*4]			; ds[(fi+1)*N+i...i+4]
		subps xmm2, xmm4					; x-muu
		mulps xmm2, xmm2 					; (x-muu)*(x-muu)
		addps xmm5, xmm2					; ximuu+=(x-muu)*(x-muu) (fi+1)

		add eax, 4							; i+=4
		add edx, 4
		add esi, 4
		jmp forN_rcf

resN_rcf:
		mov ebx, [ebp+12]					; N%p
		cmp ebx, 0
		je endN_rcf
		neg ebx								; -N%p
		add ebx, 4							; p-N%p
		sub edx, ebx						; fi*N+i -  (p-N%p)
		sub esi, ebx						; (fi+1)*N+i -  (p-N%p)
		movups xmm2, [ecx+edx*4]			; ds[fi*N+i...i+4]
		subps xmm2, xmm0					; x-muu (fi)

		movups xmm3, [ecx+esi*4]			; ds[(fi+1)*N+i...i+4]
		subps xmm3, xmm4					; x-muu (fi+1)

		movups xmm0, [maskres]				; 1 1 1 0

		cmp ebx, 2							; p-N%p == 2
		jne resN1_rcf
		insertps xmm0, xmm0, 0x3			; 1 1 0 0

resN1_rcf:
		cmp ebx, 3							; p-N%p == 3
		jne calc_resN_rcf
		insertps xmm0, xmm0, 0x7			; 1 0 0 0

calc_resN_rcf:
		andps xmm2, xmm0					; x-muu masked (fi)
		andps xmm3, xmm0					; x-muu masked (fi+1)
		mulps xmm2, xmm2 					; (x-muu)*(x-muu) (fi)
		addps xmm1, xmm2					; ximuu+=(x-muu)*(x-muu) (fi)
		mulps xmm3, xmm3 					; (x-muu)*(x-muu) (fi+1)
		addps xmm5, xmm3					; ximuu+=(x-muu)*(x-muu) (fi+1)

endN_rcf:	
		haddps xmm1, xmm1
		haddps xmm1, xmm1
		sqrtss xmm1, xmm1					; sqrt(ximuu) (fi)
		movups xmm2, [maskfloat]			; 0 0 0 1.0
		divss xmm2, xmm1					; 0 0 0 1.0/sqrt(ximuu) (fi)

		haddps xmm5, xmm5
		haddps xmm5, xmm5
		sqrtss xmm5, xmm5					; sqrt(ximuu) (fi+1)
		movups xmm3, [maskfloat]			; 0 0 0 1.0
		divss xmm3, xmm5					; 0 0 0 1.0/sqrt(ximuu) (fi+1)

		mov edx, [ebp+24]					; strcf*
		mulss xmm2, [edx + edi*4]
		movss [edx + edi*4], xmm2			; strcf[fi]
		mulss xmm3, [edx + edi*4 + 4]
		movss [edx + edi*4 + 4], xmm3		; strcf[fi+1]

		comiss xmm2, xmm3					; strcf[fi] > strcf[fi+1]
		jbe cmp_rcf
		comiss xmm2, xmm7					; strcf[fi] > rcfMax
		jbe endD_rcf
		pop esi
		push edi							; fmax = fi
		movss xmm7, xmm2					; rcfMax = strcf[fi]
		jmp endD_rcf

cmp_rcf:
		comiss xmm3, xmm7					; strcf[fi+1] > rcfMax
		jbe endD_rcf
		inc edi
		pop esi
		push edi							; fmax = fi+1
		dec edi
		movss xmm7, xmm3					; rcfMax = strcf[fi+1]

endD_rcf:
		mov esi, [ebp + 8]					; d%p
		neg esi 							; -d%p
		mov ecx, [ebp + 20]					; input
		mov ecx, [ecx + 24]					; d
		add esi, ecx						; d-d%p
		add edi, 2							; fi+2
		jmp for_rcf

resD_rcf:
		mov ebx, [ebp+8]					; d%p
		cmp ebx, 0
		je end_rcf
		mov edx, [ebp+16]					; mu*
		mov ecx, [ebp+20]					; input
		mov ecx, [ecx+24]					; d
		shl ecx, 1							; 2*d
		add ecx, edi						; 2*d+fi
		movss xmm0, [edx+ecx*4]				; x x x muu (fi)
		shufps xmm0, xmm0, 0				; muu[] (fi)

		xorps xmm1, xmm1					; ximuu (fi)
		xor eax, eax						; contatore for innestato
		mov ecx, [ebp+20]					; input
		mov edx, [ecx+20]					; N
		mov ebx, [ebp+12]					; N%p
		neg ebx								; -N%p
		add ebx, edx						; N-N%p
		mov ecx, [ecx]						; ds
		imul edx, edi						; fi*N

resD_forN_rcf: 
		cmp eax, ebx						; i < N-N%p
		jge resD_resN_rcf
		movups xmm2, [ecx+edx*4]			; ds[fi*N+i...i+4]
		subps xmm2, xmm0					; x-muu
		mulps xmm2, xmm2 					; (x-muu)*(x-muu)
		addps xmm1, xmm2					; ximuu+=(x-muu)*(x-muu) (fi)
		add eax, 4							; i+=4
		add edx, 4
		jmp resD_forN_rcf

resD_resN_rcf:
		mov ebx, [ebp+12]					; N%p
		cmp ebx, 0
		je resD_endN_rcf
		neg ebx								; -N%p
		add ebx, 4							; p-N%p
		sub edx, ebx						; fi*N+i -  (p-N%p)
		movups xmm2, [ecx+edx*4]			; ds[fi*N+i...i+4]
		subps xmm2, xmm0					; x-muu (fi)
		movups xmm0, [maskres]				; 1 1 1 0
		cmp ebx, 2							; p-N%p == 2
		jne resD_resN1_rcf
		insertps xmm0, xmm0, 0x3			; 1 1 0 0

resD_resN1_rcf:
		cmp ebx, 3							; p-N%p == 3
		jne resD_calc_resN_rcf
		insertps xmm0, xmm0, 0x7			; 1 0 0 0

resD_calc_resN_rcf:
		andps xmm2, xmm0					; x-muu masked (fi)

		mulps xmm2, xmm2 					; (x-muu)*(x-muu) (fi)
		addps xmm1, xmm2					; ximuu+=(x-muu)*(x-muu) (fi)

resD_endN_rcf:	
		haddps xmm1, xmm1
		haddps xmm1, xmm1
		sqrtss xmm1, xmm1					; sqrt(ximuu) (fi)
		movups xmm2, [maskfloat]			; 0 0 0 1.0
		divss xmm2, xmm1					; 0 0 0 1.0/sqrt(ximuu) (fi)
		mov edx, [ebp+24]					; strcf*
		mulss xmm2, [edx + edi*4]
		movss [edx + edi*4], xmm2			; strcf[fi]
		comiss xmm2, xmm7					; strcf[fi] > rcfMax
		jbe end_rcf
		pop esi
		push edi							; fmax = fi
		movss xmm7, xmm2					; rcfMax = strcf[fi]

end_rcf:	
		pop eax								; return fmax 
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante

global rffP4

rffP4:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; rffP4(fx, fy, n%p, mu*, input*)
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		
		; ------------------------------------------------------------
		; Elaborazione
		; ------------------------------------------------------------
		mov ecx, [ebp+24]					; input
		mov ecx, [ecx+20]					; N
		mov esi, [ebp+16]					; N%p
		neg esi 							; -N%p
		add esi, ecx						; N-N%p
		mov eax, [ebp+8] 					; fx
		imul eax, ecx						; fx*N
		mov ebx, [ebp+12]					; fy
		imul ebx, ecx						; fy*N
		mov edx, [ebp + 20]					; mu*
		mov ecx, [ebp + 24]					; input*
		mov ecx, [ecx + 24]					; input->d 
		shl ecx, 1							; 2*d
		mov edi, ecx						; copy 2*d
		add ecx, [ebp + 8]					; 2*d+fx
		movss xmm3, [edx+ecx*4]				; x x x mux
		shufps xmm3, xmm3, 0				; mux[]
		add edi, [ebp + 12]					; 2*d+fy
		movss xmm4, [edx+edi*4] 			; x x x muy
		shufps xmm4, xmm4, 0		 		; muy[]
		mov edx, [ebp+24]					; input
		mov edx, [edx]						; ds
		xor edi, edi
		xorps xmm5, xmm5
		xorps xmm6, xmm6
		xorps xmm7, xmm7

for_rff:
		cmp edi, esi						; i < N-N%p 
		jge res_rff

		movups xmm0, [edx+eax*4]			; ds[fx*N+i...i+4]
		subps xmm0, xmm3					; x-mux
		
		movups xmm1, [edx+ebx*4]			; ds[fy*N+i...i+4]
		subps xmm1, xmm4					; y-muy

		movups xmm2, xmm0 					; copy x-mux
		mulps xmm2, xmm1 					; (x-mux)*(y-mux)
		addps xmm5, xmm2					; num num num num

		mulps xmm0, xmm0 					; (x-mux)^2
		mulps xmm1, xmm1 					; (y-muy)^2
		addps xmm6, xmm0					; den1 den1 den1 den1 
		addps xmm7, xmm1					; den2 den2 den2 den2

		add edi, 4
		add eax, 4 							; fx*N + i
		add ebx, 4 							; fy*N + i
		jmp for_rff

res_rff:
		mov esi, [ebp+16]					; N%p
		cmp esi, 0
		je end_rff
		neg esi 							; -N%p
		add esi, 4							; p-N%p
		sub eax, esi 						; fx*N+i -  (p-N%p)
		sub ebx, esi 						; fy*N+i -  (p-N%p)
		
		movups xmm0, [edx+eax*4]			; ds[fx*N+i...i+4]
		subps xmm0, xmm3					; x-mux
		movups xmm1, [edx+ebx*4]			; ds[fy*N+i...i+4]
		subps xmm1, xmm4					; y-muy

		movups xmm3, [maskres]				; 1 1 1 0
		cmp esi, 2							; p-N%p == 2
		jne res1_rff
		insertps xmm3, xmm3, 0x3			; 1 1 0 0
		jmp calc_res_rff

res1_rff:	
		cmp esi, 3							; p-N%p == 3
		jne calc_res_rff
		insertps xmm3, xmm3, 0x7			; 1 0 0 0

calc_res_rff:
		andps xmm0, xmm3					; x-mux masked
		andps xmm1, xmm3					; y-muy masked

		movups xmm2, xmm0 					; copy x-mux
		mulps xmm2, xmm1 					; (x-mux)*(y-mux)
		addps xmm5, xmm2					; num num num num

		mulps xmm0, xmm0 					; (x-mux)^2
		mulps xmm1, xmm1 					; (y-muy)^2
		addps xmm6, xmm0					; den1 den1 den1 den1 
		addps xmm7, xmm1					; den2 den2 den2 den2

end_rff: 	
		haddps xmm5, xmm5
		haddps xmm5, xmm5					; num
		haddps xmm6, xmm6
		haddps xmm6, xmm6					; den1
		haddps xmm7, xmm7
		haddps xmm7, xmm7					; den2
		mulss xmm6, xmm7					; den1*den2

		movups xmm3, [maskabs]				; load mask abs
		andps xmm5, xmm3					; |num|
		sqrtss xmm6, xmm6					; sqrt(den1*den2)
		divss xmm5, xmm6					; rff

		movss [mem_rff], xmm5				; store rff
		emms								; switch to x87
		fld dword [mem_rff]					; return rff
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante