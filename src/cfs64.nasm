; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
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
;     nasm -f elf64 regression64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati
	maskres3 	dq 	0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF
	maskres2 	dq 	0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF
	maskres1 	dq 	0, 0, 0, 0xFFFFFFFFFFFFFFFF

	maskfloat 	dq 	0x3ff0000000000000, 0, 0, 0	
	maskabs		dq 	0x7FFFFFFFFFFFFFFF, 0, 0, 0		

section .bss			; Sezione contenente dati non inizializzati

alignb 32
sc		resq		1

section .text			; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
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
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

; ------------------------------------------------------------
; Funzione prova
; ------------------------------------------------------------
global prova

msg	db 'sc:',32,0
nl	db 10,0

prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		; ------------------------------------------------------------
		; I parametri sono passati nei registri
		; ------------------------------------------------------------
		; rdi = indirizzo della struct input
		
		; esempio: stampa input->sc
        ; [RDI] input->ds; 			// dataset
		; [RDI + 8] input->labels; 	// etichette
		; [RDI + 16] input->out;	// vettore contenente risultato dim=(k+1)
		; [RDI + 24] input->sc;		// score dell'insieme di features risultato
		; [RDI + 32] input->k; 		// numero di features da estrarre
		; [RDI + 36] input->N;		// numero di righe del dataset
		; [RDI + 40] input->d;		// numero di colonne/feature del dataset
		; [RDI + 44] input->display;
		; [RDI + 48] input->silent;
		;VMOVSD		ymm0, [RDI+24]
		;VMOVSD		[sc], ymm0
		;prints 		msg
		;printsd		sc
		;prints 		nl
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante

global mediaP4

mediaP4:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali
		
		; ------------------------------------------------------------
		; Elaborazione
		; ------------------------------------------------------------
		xor rax, rax
		xor rbx, rbx
		mov eax, [rcx + 36]			; input->N
		mov ebx, esi				; N%p 
		neg ebx						; -N%p 
		add ebx, eax				; N - N%p
		imul eax, edi				; N * fx
		mov r10, [rcx + 8]			; input->labels*
		mov r11, [rcx]				; input->ds*
		xor r15, r15				; i=0
		vxorpd ymm5, ymm5			; n1
		vxorpd ymm6, ymm6			; m1
		vxorpd ymm7, ymm7			; m0
		vxorpd ymm3, ymm3			; 0 0 0 0

for_media:
		cmp r15, rbx				; i < N-N%p
		jnl	res_media

		vmovupd ymm0, [r11 + rax*8]	; ds[fx*N+i,...,fx*N+i+3]
		vmovupd ymm1, ymm0			; copy ds[fx*N+i,...,fx*N+i+3]
		vmovupd ymm2, [r10 + r15*8]	; labels[i,...,i+3]

		vaddpd ymm5, ymm2			; accumula labels[i]=1
		vcmpneqpd ymm2, ymm3		; cmp not eq con mask di 0
		
		vandpd ymm0, ymm2			; ds and labels (c=1)
		vandnpd ymm2, ymm1			; ds and not labels (c=0)
		
		vaddpd ymm6, ymm0			; m1 m1 m1 m1 
		vaddpd ymm7, ymm2			; m0 m0 m0 m0

		add r15, 4					; i+=p
		add rax, 4					; (fx*N+i)+4
		jmp for_media

res_media:
		mov r14, rsi				; N%p 
		cmp r14, 0
		je end_media
		neg r14 					; -N%p
		add r14, 4					; p-N%p

		sub rax, r14 				; fx*N+i -  (p-N%p)
		sub r15, r14				; labels[i] -  (p-N%p)

		vmovupd ymm0, [r11 + rax*8]		; ds[fx*N+i...i+3]
		vmovupd ymm1, ymm0				; copy ds[fx*N+i,...,fx*N+i+3]
		vmovupd ymm2, [r10 + r15*8]		; labels[i,...,i+3]

		vmovupd ymm3, [maskres3]		; 1 1 1 0
		
		cmp r14, 2					 ; p-N%p == 2
		jne res1_media
		vmovupd ymm3, [maskres2] 	 ; 1 1 0 0 11|11|00|00
		jmp calc_res_media
		
res1_media:	
		cmp r14, 3					; p-N%p == 3
		jne calc_res_media
		vmovupd ymm3, [maskres1]	; 1 0 0 0 11|00|00|00

calc_res_media:
		vandpd ymm0, ymm3			; masked ds[fx*N+i...i+4]
		vandpd ymm1, ymm3			; masked ds[fx*N+i...i+4]
		vandpd ymm2, ymm3			; masked n1
		vaddpd ymm5, ymm2			; accumula labels[i]=1
		vxorpd ymm4, ymm4			; 0 0 0 0
		vcmpneqpd ymm2, ymm4		; cmp not eq con mask di 0
		vandpd ymm2, ymm3			; masked not labels

		vandpd ymm0, ymm2			; ds and labels (c=1)
		vandnpd ymm2, ymm1			; ds and not labels (c=0)
		
		vaddpd ymm6, ymm0			; m1 m1 m1 m1
		vaddpd ymm7, ymm2			; m0 m0 m0 m0

end_media:	
		vxorpd ymm0, ymm0
		vxorpd ymm1, ymm1

		vhaddpd ymm5, ymm5			 ; x n1 x n1
		vextractf128 xmm3, ymm5, 1			
		vaddpd ymm5, ymm3			; x x x n1

		vhaddpd ymm6, ymm6			; x m1 x m1
		vextractf128 xmm3, ymm6, 1			
		vaddpd ymm6, ymm3			; x x x m1

		vhaddpd ymm7, ymm7			; x m0 x m0
		vextractf128 xmm3, ymm7, 1			
		vaddpd ymm7, ymm3			; x x x m0

		vmovupd ymm0, ymm6			; copy m1
		vdivpd ymm0, ymm5			; x x x m1/n1=mu1

		cvtsi2sd xmm1, [rcx + 36]	; input->N Convert_Integer_To_Double_Precision_Floating_Point
		
		vmovupd ymm4, ymm1			; copy N

		vaddpd ymm6, ymm7			; x x x m1+m0
		vdivpd ymm6, ymm1			; x x x (m1+m0)/N=mu

		vsubpd ymm1, ymm5 			; x x x N-n1(=n0)
		vdivpd ymm7, ymm1			; x x x m0/n0=mu0

		xor rax, rax
		mov eax, [rcx + 40]			; input->d

		movsd [rdx + rdi*8], xmm7, 	; mu[fx]=mu0
		add rax, rdi				; fx+d
		movsd [rdx + rax*8], xmm0	; mu[fx+d]=mu1
		sub rax, rdi				; d
		add rax, rax				; 2d
		add rax, rdi				; fx+2*d
		movsd [rdx + rax*8], xmm6	; mu[fx+2*d]=mu

		vsubpd ymm7, ymm0			; mu0-mu1
		vmovupd ymm3, [maskabs]		; load mask abs
		vandpd ymm7, ymm3			; |mu0-mu1|
		vdivpd ymm7, ymm4			; |mu0-mu1|/N

		vsqrtpd ymm1, ymm1			; sqrt(n0)
		vsqrtpd ymm5, ymm5			; sqrt(n1)

		vsubpd ymm4, [maskfloat]	; x x x N-1
		vsqrtpd ymm4, ymm4			; sqrt(N-1)

		vmulpd ymm7, ymm1			; |mu0-mu1|/N * sqrt(n0)
		vmulpd ymm7, ymm5			; |mu0-mu1|/N * sqrt(n0) * sqrt(n1)
		vmulpd ymm7, ymm4			; |mu0-mu1|/N * sqrt(n0) * sqrt(n1) * sqrt(N-1)

		movsd [r8 + rdi*8], xmm7	; strcf[fx] = |mu0-mu1|/N * sqrt(n0) * sqrt(n1) * sqrt(N-1)
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante


global rcfP4

rcfP4:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		push rbx
		push rcx
		push rdx
		push rsi
		push rdi
		push r8
		push r9
		push r10
		push r11
		push r12
		push r13
		push r14
		push r15
		sub    rsp,0x20
  		vmovups [rsp],ymm0
  		sub    rsp,0x20
   		vmovups [rsp],ymm1
   		sub    rsp,0x20
   		vmovups [rsp],ymm2
		sub    rsp,0x20
   		vmovups [rsp],ymm3
   		sub    rsp,0x20
		vmovups [rsp],ymm4
		sub    rsp,0x20
		vmovups [rsp],ymm5
		sub    rsp,0x20
		vmovups [rsp],ymm6
		sub    rsp,0x20
		vmovups [rsp],ymm7
		sub    rsp,0x20
		vmovups [rsp],ymm8
		sub    rsp,0x20
		vmovups [rsp],ymm9
		sub    rsp,0x20
		vmovups [rsp],ymm10
		sub    rsp,0x20
		vmovups [rsp],ymm11
		sub    rsp,0x20
		vmovups [rsp],ymm12
		sub    rsp,0x20
		vmovups [rsp],ymm13
		sub    rsp,0x20
		vmovups [rsp],ymm14
		sub    rsp,0x20
		vmovups [rsp],ymm15
		; ------------------------------------------------------------
		; Elaborazione
		; ------------------------------------------------------------
		xor rax, rax
		xor rbx, rbx
		mov eax, [rcx + 40]					; d
		mov r9, rdi							; d%p
		neg r9 								; -d%p
		add r9, rax							; d-d%p
		vxorpd ymm7, ymm7					; rcfMax
		push 0								; fmax
		xor r10, r10						; fi=0
		mov r14, [rcx]						; ds*

for_rcf: 
		cmp r10, r9							; fi < d-d%p 
		jge resD_rcf
		xor rax, rax
		mov eax, [rcx + 40]					; d
		shl rax, 1							; 2*d
		add rax, r10						; 2*d+fi
		vbroadcastsd ymm0, [rdx + rax*8]		; muu muu muu muu (fi)
		vbroadcastsd ymm4, [rdx + rax*8 + 8]	; muu muu muu muu (fi+1)

		vxorpd ymm1, ymm1					; ximuu (fi)
		vxorpd ymm5, ymm5					; ximuu (fi+1)

		xor r11, r11						; i=0 
		mov ebx, [rcx + 36]					; N
		mov r12, rsi						; N%p
		neg r12								; -N%p
		add r12, rbx						; N-N%p

		mov r13, r10						; copy fi
		inc r13								; fi+1
		imul r13, rbx						; (fi+1)*N
		imul rbx, r10						; fi*N

forN_rcf: 
		cmp r11, r12						; i < N-N%p
		jge resN_rcf

		vmovupd ymm2, [r14 + rbx*8]			; ds[fi*N+i...i+4]
		vsubpd ymm2, ymm0					; x-muu
		
		vmulpd ymm2, ymm2 					; (x-muu)*(x-muu)
		vaddpd ymm1, ymm2					; ximuu+=(x-muu)*(x-muu) (fi)

		vmovupd ymm2, [r14 + r13*8]			; ds[(fi+1)*N+i...i+4]
		vsubpd ymm2, ymm4					; x-muu
		
		vmulpd ymm2, ymm2 					; (x-muu)*(x-muu)
		vaddpd ymm5, ymm2					; ximuu+=(x-muu)*(x-muu) (fi+1)

		add r11, 4							; i+=4
		add rbx, 4							; fi*N + 4
		add r13, 4							; (fi+1)*N + 4
		jmp forN_rcf

resN_rcf:
		cmp rsi, 0
		je endN_rcf
		mov r12, rsi						; N%p
		neg r12								; -N%p
		add r12, 4							; p-N%p

		sub rbx, r12						; fi*N+i -  (p-N%p)
		sub r13, r12						; (fi+1)*N+i -  (p-N%p)

		vmovupd ymm2, [r14 + rbx*8]			; ds[fi*N+i...i+4]
		vsubpd ymm2, ymm0					; x-muu (fi)

		vmovupd ymm3, [r14 + r13*8]			; ds[(fi+1)*N+i...i+4]
		vsubpd ymm3, ymm4					; x-muu (fi+1)

		vmovupd ymm0, [maskres3]			; 1 1 1 0

		cmp r12, 2							; p-N%p == 2
		jne resN1_rcf
		vmovupd ymm0, [maskres2]			; 1 1 0 0 11|11|00|00

resN1_rcf:
		cmp r12, 3							; p-N%p == 3
		jne calc_resN_rcf
		vmovupd ymm0, [maskres1]			; 1 0 0 0 11|00|00|00

calc_resN_rcf:
		vandpd ymm2, ymm0					; x-muu masked (fi)
		vandpd ymm3, ymm0					; x-muu masked (fi+1)

		vmulpd ymm2, ymm2 					; (x-muu)*(x-muu) (fi)
		vaddpd ymm1, ymm2					; ximuu+=(x-muu)*(x-muu) (fi)

		vmulpd ymm3, ymm3 					; (x-muu)*(x-muu) (fi+1)
		vaddpd ymm5, ymm3					; ximuu+=(x-muu)*(x-muu) (fi+1)

endN_rcf:	
		vhaddpd ymm1, ymm1
		vextractf128 xmm3, ymm1, 1			
		vaddpd ymm1, ymm3					; x x x res
		vsqrtpd ymm1, ymm1					; sqrt(ximuu) (fi)
		vmovupd ymm2, [maskfloat]			; 0 0 0 1.0
		vdivpd ymm2, ymm1					; 0 0 0 1.0/sqrt(ximuu) (fi)

		vhaddpd ymm5, ymm5
		vextractf128 xmm3, ymm5, 1			
		vaddpd ymm5, ymm3					; x x x res
		vsqrtpd ymm5, ymm5					; sqrt(ximuu) (fi+1)
		
		vmovupd ymm3, [maskfloat]			; 0 0 0 1.0
		vdivpd ymm3, ymm5					; 0 0 0 1.0/sqrt(ximuu) (fi+1)

		mulsd xmm2, [r8 + r10*8]
		movsd [r8 + r10*8], xmm2			; strcf[fi]
		mulsd xmm3, [r8 + r10*8 + 8]
		movsd [r8 + r10*8 + 8], xmm3		; strcf[fi+1]

		comisd xmm2, xmm3					; strcf[fi] > strcf[fi+1]
		jbe cmp_rcf
		comisd xmm2, xmm7					; strcf[fi] > rcfMax
		jbe endD_rcf
		pop r15								
		push r10							; fmax = fi
		vmovupd ymm7, ymm2					; rcfMax = strcf[fi]
		jmp endD_rcf

cmp_rcf:
		comisd xmm3, xmm7					; strcf[fi+1] > rcfMax
		jbe endD_rcf
		inc r10
		pop r15
		push r10							; fmax = fi+1
		dec r10
		vmovupd ymm7, ymm3					; rcfMax = strcf[fi+1]

endD_rcf:
		add r10, 2							; fi+2
		jmp for_rcf

resD_rcf:
		cmp rdi, 0
		je end_rcf
		mov r9, rdi							; d%p
		xor rax, rax
		mov eax, [rcx + 40]					; d
		shl rax, 1							; 2*d
		add rax, r10						; 2*d+fi

		vbroadcastsd ymm0, [rdx + rax*8]	; muu muu muu muu (fi)

		vxorpd ymm1, ymm1					; ximuu (fi)
		xor r11, r11						; i=0
		xor rbx, rbx
		mov ebx, [rcx + 36]					; N
		mov r12, rsi						; N%p
		neg r12								; -N%p
		add r12, rbx						; N-N%p
		imul rbx, r10						; fi*N

resD_forN_rcf: 
		cmp r11, r12						; i < N-N%p
		jge resD_resN_rcf

		vmovupd ymm2, [r14 + rbx*8]			; ds[fi*N+i...i+4]
		vsubpd ymm2, ymm0					; x-muu
		vmulpd ymm2, ymm2 					; (x-muu)*(x-muu)
		vaddpd ymm1, ymm2					; ximuu+=(x-muu)*(x-muu) (fi)

		add r11, 4							; i+=4
		add rbx, 4							; fi*N + 4
		jmp resD_forN_rcf

resD_resN_rcf:
		cmp rsi, 0
		je resD_endN_rcf
		xor rax, rax
		mov rax, rsi						; N%p
		neg rax								; -N%p
		add rax, 4							; p-N%p

		sub rbx, rax						; fi*N+i -  (p-N%p)

		vmovupd ymm2, [r14 + rbx*8]			; ds[fi*N+i...i+4]
		vsubpd ymm2, ymm0					; x-muu (fi)

		vmovupd ymm0, [maskres3]			; 1 1 1 0

		cmp rax, 2							; p-N%p == 2
		jne resD_resN1_rcf
		vmovupd ymm0, [maskres2] 			; 1 1 0 0 11|11|00|00

resD_resN1_rcf:
		cmp rax, 3							; p-N%p == 3
		jne resD_calc_resN_rcf
		vmovupd ymm0, [maskres1]			; 1 0 0 0 11|00|00|00

resD_calc_resN_rcf:
		vandpd ymm2, ymm0					; x-muu masked (fi)

		vmulpd ymm2, ymm2 					; (x-muu)*(x-muu) (fi)
		vaddpd ymm1, ymm2					; ximuu+=(x-muu)*(x-muu) (fi)

resD_endN_rcf:	
		vhaddpd ymm1, ymm1
		vextractf128 xmm3, ymm1, 1			
		vaddpd ymm1, ymm3					; x x x res
		vsqrtpd ymm1, ymm1					; sqrt(ximuu) (fi)
		vmovupd ymm2, [maskfloat]			; 0 0 0 1.0
		vdivpd xmm2, xmm1					; 0 0 0 1.0/sqrt(ximuu) (fi)

		mulsd xmm2, [r8 + r10*8]
		movsd [r8 + r10*8], xmm2			; strcf[fi]

		comisd xmm2, xmm7					; strcf[fi] > rcfMax
		jbe end_rcf
		pop r15
		push r10							; fmax = fi
		vmovupd ymm7, ymm2					; rcfMax = strcf[fi]

end_rcf:	
		pop rax								; return fmax 
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		vmovups ymm15,[rsp]
   		add    rsp,0x20
   		vmovups ymm14,[rsp]
   		add    rsp,0x20
		vmovups ymm13,[rsp]
		add    rsp,0x20
		vmovups ymm12,[rsp]
		add    rsp,0x20
		vmovups ymm11,[rsp]
		add    rsp,0x20
		vmovups ymm10,[rsp]
		add    rsp,0x20
		vmovups ymm9,[rsp]
		add    rsp,0x20
		vmovups ymm8,[rsp]
		add    rsp,0x20
		vmovups ymm7,[rsp]
		add    rsp,0x20
		vmovups ymm6,[rsp]
		add    rsp,0x20
		vmovups ymm5,[rsp]
		add    rsp,0x20
		vmovups ymm4,[rsp]
		add    rsp,0x20
		vmovups ymm3,[rsp]
		add    rsp,0x20
		vmovups ymm2,[rsp]
		add    rsp,0x20
		vmovups ymm1,[rsp]
		add    rsp,0x20
		vmovups ymm0,[rsp]
		add    rsp,0x20
		pop    r15
		pop    r14
		pop    r13
		pop    r12
		pop    r11
		pop    r10
		pop    r9
		pop    r8
		pop    rdi
		pop    rsi
		pop    rdx
		pop    rcx
		pop    rbx
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp			; ripristina il Base Pointer
		ret					; torna alla funzione C chiamante
		
global rffP4

rffP4:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		push rax
		push rbx
		push rcx
		push rdx
		push rsi
		push rdi
		push r8
		push r9
		push r10
		push r11
		push r12
		push r13
		push r14
		push r15
  		sub    rsp,0x20
   		vmovups [rsp],ymm1
   		sub    rsp,0x20
   		vmovups [rsp],ymm2
		sub    rsp,0x20
   		vmovups [rsp],ymm3
   		sub    rsp,0x20
		vmovups [rsp],ymm4
		sub    rsp,0x20
		vmovups [rsp],ymm5
		sub    rsp,0x20
		vmovups [rsp],ymm6
		sub    rsp,0x20
		vmovups [rsp],ymm7
		sub    rsp,0x20
		vmovups [rsp],ymm8
		sub    rsp,0x20
		vmovups [rsp],ymm9
		sub    rsp,0x20
		vmovups [rsp],ymm10
		sub    rsp,0x20
		vmovups [rsp],ymm11
		sub    rsp,0x20
		vmovups [rsp],ymm12
		sub    rsp,0x20
		vmovups [rsp],ymm13
		sub    rsp,0x20
		vmovups [rsp],ymm14
		sub    rsp,0x20
		vmovups [rsp],ymm15
		; ------------------------------------------------------------
		; Elaborazione
		; ------------------------------------------------------------
		xor rax, rax
		xor rbx, rbx
		mov eax, [r8 + 36]					; N
		mov r9, rdi 						; fx
		imul r9, rax						; fx*N
		mov r10, rsi						; fy
		imul r10, rax						; fy*N
		mov ebx, edx						; N%p
		neg ebx 							; -N%p
		add ebx, eax						; N-N%p
		xor rax, rax
		mov eax, [r8 + 40]					; input->d 
		shl eax, 1							; 2*d
		mov r12, rax						; copy 2*d 
		add rax, rdi						; 2*d+fx
		vbroadcastsd ymm3, [rcx + rax*8]	; mux mux mux mux
		add r12, rsi						; 2*d+fy
		vbroadcastsd ymm4, [rcx + r12*8] 	; muy muy muy muy
		mov r13, [r8]
		xor r11, r11
		vxorpd ymm5, ymm5
		vxorpd ymm6, ymm6
		vxorpd ymm7, ymm7

for_rff:
		cmp r11, rbx						; i < N-N%p 
		jge res_rff

		vmovupd ymm0, [r13 + r9*8]			; ds[fx*N+i...i+4]
		vsubpd ymm0, ymm3					; x-mux
		
		vmovupd ymm1, [r13 + r10*8]			; ds[fy*N+i...i+4]
		vsubpd ymm1, ymm4					; y-muy

		vmovupd ymm2, ymm0 					; copy x-mux
		vmulpd ymm2, ymm1 					; (x-mux)*(y-mux)
		vaddpd ymm5, ymm2					; num num num num

		vmulpd ymm0, ymm0 					; (x-mux)^2
		vmulpd ymm1, ymm1 					; (y-muy)^2
		vaddpd ymm6, ymm0					; den1 den1 den1 den1 
		vaddpd ymm7, ymm1					; den2 den2 den2 den2

		add r11, 4
		add r9, 4 							; fx*N + i
		add r10, 4 							; fy*N + i
		jmp for_rff

res_rff:
		mov ebx, edx						; N%p
		cmp rbx, 0
		je end_rff
		neg rbx 							; -N%p
		add rbx, 4							; p-N%p

		sub r9, rbx 						; fx*N+i -  (p-N%p)
		sub r10, rbx 						; fy*N+i -  (p-N%p)

		vmovupd ymm0, [r13 + r9*8]			; ds[fx*N+i...i+4]
		vsubpd ymm0, ymm3					; x-mux
		vmovupd ymm1, [r13 + r10*8]			; ds[fy*N+i...i+4]
		vsubpd ymm1, ymm4					; y-muy

		vmovupd ymm3, [maskres3]			; 1 1 1 0
		
		cmp rbx, 2					 		; p-N%p == 2
		jne res1_rff
		vmovupd ymm3, [maskres2]	 		; 1 1 0 0 11|11|00|00
		jmp calc_res_rff

res1_rff:	
		cmp rbx, 3							; p-N%p == 3
		jne calc_res_rff
		vmovupd ymm3, [maskres1]			; 1 0 0 0 11|00|00|00

calc_res_rff:
		vandpd ymm0, ymm3					; x-mux masked
		vandpd ymm1, ymm3					; y-muy masked

		vmovupd ymm2, ymm0 					; copy x-mux
		vmulpd ymm2, ymm1 					; (x-mux)*(y-mux)
		vaddpd ymm5, ymm2					; num num num num

		vmulpd ymm0, ymm0 					; (x-mux)^2
		vmulpd ymm1, ymm1 					; (y-muy)^2
		vaddpd ymm6, ymm0					; den1 den1 den1 den1 
		vaddpd ymm7, ymm1					; den2 den2 den2 den2

end_rff: 	
		vhaddpd ymm5, ymm5					; x num x num
		vextractf128 xmm3, ymm5, 1			
		vaddpd ymm5, ymm3					; x x x num

		vhaddpd ymm6, ymm6					; x den x den
		vextractf128 xmm3, ymm6, 1			
		vaddpd ymm6, ymm3					; x x x den

		vhaddpd ymm7, ymm7					; x den x den
		vextractf128 xmm3, ymm7, 1			
		vaddpd ymm7, ymm3					; x x x den
		vmulpd ymm6, ymm7					; x x x den1*den2

		vmovupd ymm3, [maskabs]				; load mask abs
		vandpd ymm5, ymm3					; |num|
		vsqrtpd ymm6, ymm6					; sqrt(den1*den2)
		vdivpd ymm5, ymm6					; rff

		vmovupd ymm0, ymm5					; return rff
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		vmovups ymm15,[rsp]
   		add    rsp,0x20
   		vmovups ymm14,[rsp]
   		add    rsp,0x20
		vmovups ymm13,[rsp]
		add    rsp,0x20
		vmovups ymm12,[rsp]
		add    rsp,0x20
		vmovups ymm11,[rsp]
		add    rsp,0x20
		vmovups ymm10,[rsp]
		add    rsp,0x20
		vmovups ymm9,[rsp]
		add    rsp,0x20
		vmovups ymm8,[rsp]
		add    rsp,0x20
		vmovups ymm7,[rsp]
		add    rsp,0x20
		vmovups ymm6,[rsp]
		add    rsp,0x20
		vmovups ymm5,[rsp]
		add    rsp,0x20
		vmovups ymm4,[rsp]
		add    rsp,0x20
		vmovups ymm3,[rsp]
		add    rsp,0x20
		vmovups ymm2,[rsp]
		add    rsp,0x20
		vmovups ymm1,[rsp]
		add    rsp,0x20
		pop    r15
		pop    r14
		pop    r13
		pop    r12
		pop    r11
		pop    r10
		pop    r9
		pop    r8
		pop    rdi
		pop    rsi
		pop    rdx
		pop    rcx
		pop    rbx
		pop    rax
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp			; ripristina il Base Pointer
		ret					; torna alla funzione C chiamante
