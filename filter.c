#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/**
 * Gogea Mihail 
 *  
 */



int main(int argc, char * argv[]) {
    int rank;
    int nProcesses ; 

    int i,j;
    
    //Tag-uri
    int sonda = 1;
    int ecou = 2; 
    int ecou_vid = 3;
    int tag_end = 8;

    int finish = 1;

    // Detali asupra unui nod
    int nrVecini = 0;
    int parent;

    int ZERO = 0;

    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
    
    int liniiCopil[nProcesses];
    int statistica[nProcesses];

    // Matricea de filtraj + factorul la care trebuie impartit rezultatul

    int mat_filtru[3][3];
    int Factor;

    for(i = 0; i < nProcesses; i++){
        liniiCopil[i] = 0;
        statistica[i] = 0;
    }

    FILE *f;
    f = fopen(argv[1], "r");

    // Matricile pentru topologie   

    //Intializez matricile cu 0 la inceput
    int topTree[nProcesses][nProcesses];
    for(i = 0; i < nProcesses; i++){
    	for(j = 0; j < nProcesses; j++){
    		topTree[i][j] = 0;    	}
    }

    int aux_topTree[nProcesses][nProcesses];
    for(i = 0; i < nProcesses; i++){
        for(j = 0; j < nProcesses; j++){
            aux_topTree[i][j] = 0;
        }
    }

    int null_topTree[nProcesses][nProcesses];
    for(i = 0; i < nProcesses; i++){
        for(j = 0; j < nProcesses; j++){
            null_topTree[i][j] = 0;
        }
    }


    //Citesc din fisier pentru fiecare vecin al sau a unui nod.
    while(!feof(f)){

    	int nod;
    	char cr = '\0';
    	int vecin;
    	fscanf(f, "%d", &nod);

    	if(nod == rank){
    		fscanf(f, "%c", &cr);
    		fscanf(f, "%c", &cr); 

    		while(!feof(f) && cr != '\n' ){
    			fscanf(f, "%d", &vecin);

				fscanf(f, "%c", &cr);    		

                topTree[rank][vecin] = 1;	
    		}
    		break;
    	}
    	else{
    		
    		fscanf(f, "%c", &cr);
  			while(cr != '\n')
  				fscanf(f, "%c", &cr);
    		}
    }
    fclose(f);


   
    // Crearea topologiei si eliminarea buclelor.

    if(rank == ZERO){

    	// Trimit sondaje
    	for(i = 1; i < nProcesses; i++){

    		if(topTree[rank][i] == 1){
    			
    			MPI_Send(aux_topTree, (nProcesses)*(nProcesses), MPI_INT, i, sonda, MPI_COMM_WORLD);
    			nrVecini++;
    		}
    		
    	}
    }else {
    	//Se vor primi mesaje
    	MPI_Recv(aux_topTree, (nProcesses)*(nProcesses), MPI_INT, MPI_ANY_SOURCE, sonda, MPI_COMM_WORLD, &status);

    	//Parintele va deveni sursa
        if(parent != -2)
    	   parent = status.MPI_SOURCE;

    	// Trimit sondaj mai departe la ceilalti vecini

    	for(i=0; i < nProcesses; i++){
    		if( topTree[rank][i] == 1 && i != parent ){
    			nrVecini++;
    			MPI_Send(aux_topTree, (nProcesses)*(nProcesses), MPI_INT, i, sonda, MPI_COMM_WORLD);
    		}
    	}


    }



    int e = nrVecini;

    while(e > 0){

    	MPI_Recv(aux_topTree, (nProcesses)*(nProcesses), MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    	// Daca primesc sonda voi trimite un ecou vid
    	if(status.MPI_TAG == sonda){

    		topTree[rank][status.MPI_SOURCE] = 0;
    		MPI_Send(null_topTree, (nProcesses)*(nProcesses), MPI_INT, status.MPI_SOURCE, ecou_vid, MPI_COMM_WORLD);

    	}else if(status.MPI_TAG == ecou_vid){
            e--;
            topTree[rank][status.MPI_SOURCE] = 0;
        }
        else {

    		for(i = 0; i<nProcesses; i++){
    			for(j = 0; j<nProcesses; j++){

    				topTree[i][j] = topTree[i][j] | aux_topTree[i][j];	
    				 		
    			}
    		}
    		e--;
    	}    	
    }


    // Se va trimite ecou in sus spre parinti cu raspunsul fiind topologia calculata
    if(rank != 0){
    	MPI_Send(topTree, (nProcesses)*(nProcesses), MPI_INT, parent, ecou, MPI_COMM_WORLD);
    }


    // Extind topologia in jos acuma dupa ce a fost aflata in intregime de rank0
    if(rank != 0 ){
    	MPI_Recv(topTree, (nProcesses)*(nProcesses), MPI_INT, MPI_ANY_SOURCE, ecou, MPI_COMM_WORLD, &status);
    }

    // Trimit extinderea topologiei
    for(i = 0; i < nProcesses; i++){
    	if(topTree[rank][i] == 1 && i != parent ){
    		MPI_Send(topTree, (nProcesses)*(nProcesses), MPI_INT, i, ecou, MPI_COMM_WORLD);
    	}
    }

    // Calculez din nou numarul de vecini doarece am eliminat buclele.
    nrVecini = 0;
    for(i = 0; i < nProcesses; i++){
        nrVecini += topTree[rank][i];
    }

    // Terminarea topologiei.

    char filtru[128];
    char imagineCitit[128];
    char imagineScris[128];
    int tip_filtru ;

    //********************************************** De aici incepe prelucrarea imaginilor.
    MPI_Barrier(MPI_COMM_WORLD);

    // Date pentru trimiterea segmentelor
    int nrImagini ;
    int nrLinii;

    //Numarul de segmente de trimis ( catul si restul).
    int cat;
    int rest;


    if(rank == 0){

        //Deschid fisierul imagini.in
        FILE *g, *im;
        g = fopen(argv[2], "r");

        
        fscanf(g, "%d", &nrImagini);
        char buf[112];
        int inaltime;
        int latime;
        int maximumGray;



        //Trimit la prelucrat toate imaginile din fisierul de imagini.
        for(i = 0; i < nrImagini; i++){
            char c = '\0';

            fscanf(g, "%s", filtru);
            fscanf(g, "%s", imagineCitit);
            fscanf(g, "%s", imagineScris);

             //Fisierul in care voi scrie.
            FILE *out;
            out = fopen(imagineScris, "w");
        
            //Verific tipul filtrului.
            if(filtru[0] == 'b'){
                tip_filtru = 2;
            }else if(filtru[0] == 's'){
                if(filtru[1] == 'm'){
                    tip_filtru = 1;
                }else{
                    tip_filtru = 3;
                }
            }else {
                tip_filtru = 4;
            }

            // Deschi imaginea.
            im = fopen(imagineCitit, "r");
                
            /*
            !!!!!!  Mentionez ca am presupus in urma exemplelor ca un fisier pgm
                    trebuie sa aibe la inceput un tip (ex. P2) , urmand apoi
                    cu o linie de comentariu ce incepe cu # , iar apoi pe o linie
                    latimea si inaltimea matricei de pixeli , apoi numarul de maxim
                    al culori, urmand latime*inaltime pixeli.
            */

            int nr;

            // Citesc tipul imagini.
            fscanf(im,"%s", buf);
            fprintf(out, "%s\n", buf);

            //Acum sa vedem daca avem comentariu sau nu
            fscanf(im,"%s", buf);
            if(buf[0] == '#'){

                fprintf(out, "%s",buf);
                
                while(c != '\n'){
                    fscanf(im, "%c" , &c);
                   fprintf(out,"%c", c);
                }
            }else{
                fclose(im);
                im = fopen(imagineCitit, "r");
                fscanf(im,"%s", buf);
            }

            
            
            fscanf(im, "%d", &latime);
            fprintf(out, "%d ",latime );
            fscanf(im, "%d", &inaltime);
            fprintf(out, "%d\n",inaltime );   
            fscanf(im, "%d", &maximumGray);
            fprintf(out, "%d\n", maximumGray );
            
            
            int k;

            //Matricea pe care o voi trimite la copii
            int *vectorMatrice;
            vectorMatrice = (int *)calloc( (inaltime+2) * (latime+2), sizeof(int));

            int ind = 0;

            for(j = 0; j < (inaltime+2) ; j++){
                for(k = 0; k < (latime+2); k++){
                    if(j == 0 || k == 0 || (j == inaltime + 1) || (k == latime + 1) ){
                        
                        vectorMatrice[ind++] = 0;
                    }else {
                        fscanf(im, "%d" , &vectorMatrice[ind++]);
                    }
               
                }
          
            }

            fclose(im);

            // Caut numarul de copii si le asignez un numar de linii de prelucrat

            cat = inaltime / nrVecini;
            rest = inaltime % nrVecini;

            // Asignez numarul de linii ce trebuie procesat de copii
            for(j = 0 ; j < nProcesses; j++){
                if(topTree[0][j] == 1){
                    liniiCopil[j] = cat;
                }
            }

            if(cat == 0){
            	for(j = 0 ; j < nProcesses; j++){
	                if(rest > 0){

	                    if(topTree[0][j] == 1){
	                        liniiCopil[j] += 1;
	                        rest--;
	                    }
	                }else{
	                    break;
	                }
                
            	}
            }else{
            	for(j = nProcesses - 1; j >= 0; j--){
            		if(topTree[0][j] == 1){
            			liniiCopil[j] += rest;
            			break; 
            		}
            	}
            }



            //Trimit numarul de linii ce trebuie procesat copiilor
            for(j = 0 ; j < nProcesses; j++){
                if(topTree[0][j] == 1){
                    

                    //Trimit tipul de filtru ce va fii aplicat si maximumGray
                    MPI_Send(&tip_filtru, 1, MPI_INT, j, 5, MPI_COMM_WORLD);
                    MPI_Send(&maximumGray, 1, MPI_INT, j, 5, MPI_COMM_WORLD);

                    MPI_Send(&liniiCopil[j], 1, MPI_INT, j, 5, MPI_COMM_WORLD);
                    MPI_Send(&latime, 1, MPI_INT, j, 5, MPI_COMM_WORLD);

                    //Trimit de la suma liniilor de pana la el
                    int sumaLini = 0;
                    for(k = 0; k < j; k++){
                        sumaLini += liniiCopil[k];
                    }

                    MPI_Send(vectorMatrice + (latime+2)*(sumaLini), (latime+2)*(liniiCopil[j] + 2), MPI_INT, j, 5, MPI_COMM_WORLD);

                }
            }
            
            /////////////////////////////////////// Procedeu primire date ////////////////////////////////////
            
            int maximPrimire = -2;

            // Calculez maximul de linii de prelucrat la copii , ca sa creez un vector maxim.
            for(j = 0; j < nProcesses; j++){
                if(maximPrimire < liniiCopil[j])
                    maximPrimire = liniiCopil[j];
            }

            //Vector cu datele de la copii (maximul)
            int *primire;
            primire = (int *) calloc((maximPrimire)*(latime+2),sizeof(int));


            int e = nrVecini;

            //  Primesc datele de la copii
            while(e > 0){

                MPI_Recv(primire, (latime+2)*(maximPrimire), MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &status);
                e--;

                int sursa = status.MPI_SOURCE;

                int sumaLini = 0;
                    for(k = 0; k < sursa; k++){
                        sumaLini += liniiCopil[k];
                    }

                memcpy(vectorMatrice + (latime+2) + (latime+2)*(sumaLini), primire , ((latime + 2)*(liniiCopil[sursa]))*sizeof(int));

            }

            
            
            // AFisez matricea
            ind = 0;
            for(j = 0; j < (inaltime+2) ; j++){
                for(k = 0; k < (latime+2); k++){
                   
                    if(j == 0 || k == 0 || (j == inaltime + 1) || (k == latime + 1) ){
                        ind++;
                    }else{
                        fprintf(out, "%d\n", vectorMatrice[ind]);
                        ind++;
                    }
                  
                    

                    
                }
              
            }
            

            fclose(out);   
        }
        fclose(g);

        // Trimit semnalul de terminare cu tagul tag_end
        // ca copii sa stie ca s-a terminat de prelucrat imagini.
        
        for(j = 0 ; j < nProcesses; j++){

                if(topTree[0][j] == 1){
                    
                     MPI_Send(&finish, 1, MPI_INT, j, tag_end, MPI_COMM_WORLD);


                }
        }

        int e = nrVecini;

            //  Primesc statisticile de la copii
            while(e > 0){

                MPI_Recv(liniiCopil, nProcesses, MPI_INT, MPI_ANY_SOURCE, tag_end, MPI_COMM_WORLD, &status);

                for(j = 0; j < nProcesses; j++){
                    if(liniiCopil[j] != 0){
                        statistica[j] = liniiCopil[j];
                    }
                }
                e--;

            }

        //Scriu statistica in fisier statistica.out
        FILE *fStat ;
        fStat = fopen(argv[3], "w");
        for(j = 0 ; j < nProcesses; j++){

              fprintf(fStat, "%d: %d" ,j, statistica[j]);
              if(j < nProcesses - 1)
                fprintf(fStat, "\n");
        }   
        fclose(fStat); 

    }
    
 

    if(rank != 0){

        int inaltime;
        int latime;
        int maximumGray;

        // Paradigma HEARRTBEAT primesc si trimit date in continu pana cand primesc tagul de terminare.
        while(1){

             MPI_Recv(&tip_filtru, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

             if(status.MPI_TAG == tag_end){

                for(j = 0 ; j < nProcesses; j++){

                    if(topTree[rank][j] == 1 && j != parent){
                        
                         MPI_Send(&finish, 1, MPI_INT, j, tag_end, MPI_COMM_WORLD);


                    }
                }

                int e = nrVecini - 1;

                //  Primesc statisticile de la copii
                while(e > 0){

                    MPI_Recv(liniiCopil, nProcesses, MPI_INT, MPI_ANY_SOURCE, tag_end, MPI_COMM_WORLD, &status);

                    for(j = 0; j < nProcesses; j++){
                        if(liniiCopil[j] != 0){
                            statistica[j] = liniiCopil[j];
                        }
                    }
                    e--;

                }

                MPI_Send(statistica, nProcesses,  MPI_INT, parent, tag_end, MPI_COMM_WORLD);

                break;

             }else {


             
                MPI_Recv(&maximumGray, 1, MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &status);

                //Trimit numarul de imagini.
                for(j = 0 ; j < nProcesses; j++){

                        if(topTree[rank][j] == 1 && j != parent){

                             //Trimit tipul de filtru ce va fii aplicat
                            MPI_Send(&tip_filtru, 1, MPI_INT, j, 5, MPI_COMM_WORLD);
                            MPI_Send(&maximumGray, 1, MPI_INT, j, 5, MPI_COMM_WORLD);
                        }
                }

                
                

                    MPI_Recv(&inaltime, 1, MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &status);
                    MPI_Recv(&latime, 1, MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &status);


                    

                    int *matrice;
                    matrice = (int *)calloc( (inaltime+2) * (latime+2),sizeof(int));

                    MPI_Recv(matrice, (inaltime+2)*(latime+2) , MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &status);

                    if(nrVecini > 1){
                        cat = inaltime / (nrVecini-1);
                        rest = inaltime % (nrVecini-1);

                        // Asignez numarul de linii ce trebuie procesat de copii
                        for(j = 0 ; j < nProcesses; j++){
                            if(topTree[rank][j] == 1 && j != parent){
                                liniiCopil[j] = cat;
                            }
                        }

                         if(cat == 0){
			            	for(j = 0 ; j < nProcesses; j++){
				                if(rest > 0){

				                    if(topTree[rank][j] == 1 && j != parent){
				                        liniiCopil[j] += 1;
				                        rest--;
				                    }
				                }else{
				                    break;
				                }
			                
			            	}
			            }else{
			            	for(j = nProcesses - 1; j >= 0; j--){
			            		if(topTree[rank][j] == 1 && j != parent ){
			            			liniiCopil[j] += rest;
			            			break; 
			            		}
			            	}
			            }

                            //Trimit numarul de linii ce trebuie procesat copiilor
                        for(j = 0 ; j < nProcesses; j++){
                            if(topTree[rank][j] == 1 && j != parent){
                                
                                MPI_Send(&liniiCopil[j], 1, MPI_INT, j, 5, MPI_COMM_WORLD);
                                MPI_Send(&latime, 1, MPI_INT, j, 5, MPI_COMM_WORLD);

                                    //Trimit de la suma liniilor de pana la el
                                int k;
                                int sumaLini = 0;
                                for(k = 0; k < j; k++){
                                    sumaLini += liniiCopil[k];
                                }

                                MPI_Send(matrice + (latime+2)*(sumaLini), (latime+2)*(liniiCopil[j] + 2), MPI_INT, j, 5, MPI_COMM_WORLD);

                            }
                        }


                        ////////////////////////////////////// Procedeu primire ///////////////////////

                        int maximPrimire = -2;

                        for(j = 0; j < nProcesses; j++){
                            if(maximPrimire < liniiCopil[j])
                                maximPrimire = liniiCopil[j];
                        }

                        //Vector cu datele de la copii (maximul)
                        int *primire;
                        primire = (int *) calloc((maximPrimire)*(latime+2),sizeof(int));


                        int e = nrVecini - 1;

                        //  Primesc datele de la copii
                        while(e > 0){

                            MPI_Recv(primire, (latime+2)*(maximPrimire), MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &status);
                            e--;

                            int sursa = status.MPI_SOURCE;

                            int sumaLini = 0;
                            int k;
                                for(k = 0; k < sursa; k++){
                                    sumaLini += liniiCopil[k];
                                }

                            memcpy(matrice + (latime+2) + (latime+2)*(sumaLini), primire , ((latime + 2)*(liniiCopil[sursa]))*sizeof(int)) ;

                        }

                        //////////////////////////////////// Procedeu trimitere la parinte ce am primit ///////////////////////


                        MPI_Send(matrice + (latime+2), (latime+2)*(inaltime), MPI_INT, parent, 5, MPI_COMM_WORLD);


                    }else {


                        if(inaltime > 0){
                            ////////////////////// INCEP modificari asupra unui pixelii.

                            if(tip_filtru == 1){

                                mat_filtru[0][0] = 1;
                                mat_filtru[0][1] = 1;
                                mat_filtru[0][2] = 1;
                                mat_filtru[1][0] = 1;
                                mat_filtru[1][1] = 1;      
                                mat_filtru[1][2] = 1;
                                mat_filtru[2][0] = 1;
                                mat_filtru[2][1] = 1;
                                mat_filtru[2][2] = 1;
                             
                                Factor = 9;

                            }else if(tip_filtru == 2){

                                mat_filtru[0][0] = 1;  
                                mat_filtru[0][1] = 2;
                                mat_filtru[0][2] = 1;
                                mat_filtru[1][0] = 2;
                                mat_filtru[1][1] = 4;       
                                mat_filtru[1][2] = 2;
                                mat_filtru[2][0] = 1;
                                mat_filtru[2][1] = 2;
                                mat_filtru[2][2] = 1;
                               
                                Factor = 16;

                            }else if(tip_filtru == 3){
                                mat_filtru[0][0] =  0;
                                mat_filtru[0][1] =  -2;
                                mat_filtru[0][2] =  0;
                                mat_filtru[1][0] =  -2;
                                mat_filtru[1][1] =  11;        
                                mat_filtru[1][2] =  -2;
                                mat_filtru[2][0] =  0;
                                mat_filtru[2][1] =  -2;
                                mat_filtru[2][2] =  0;
                           
                                Factor = 3;

                            }else {
                                mat_filtru[0][0] =  -1;
                                mat_filtru[0][1] =  -1;
                                mat_filtru[0][2] =  -1;
                                mat_filtru[1][0] =  -1;
                                mat_filtru[1][1] =   9;     
                                mat_filtru[1][2] =  -1;
                                mat_filtru[2][0] =  -1;
                                mat_filtru[2][1] =  -1;
                                mat_filtru[2][2] =  -1;
                                
                                Factor = 1;
                            }

                            int *temp;
                            temp = (int *)calloc( (inaltime+2) * (latime+2), sizeof(int));

                            // Modific pixelul cu ajutorul matricei de filtru.
                            for(j = (latime + 2) + 1; j < (latime + 2) + (latime+2)*(inaltime) - 1; j++){ 

                                temp[j] = (( mat_filtru[0][0] * matrice[j - (latime + 2) -1] ) +
                                       ( mat_filtru[0][1] * matrice[j - (latime + 2)] ) +
                                       ( mat_filtru[0][2] * matrice[j - (latime + 2) + 1] ) +
                                       ( mat_filtru[1][0] * matrice[j - 1] ) +
                                       ( mat_filtru[1][1] * matrice[ j ] ) +             
                                       ( mat_filtru[1][2] * matrice[j + 1] ) +
                                       ( mat_filtru[2][0] * matrice[j + (latime + 2) - 1] ) +
                                       ( mat_filtru[2][1] * matrice[j + (latime + 2)] ) +
                                       ( mat_filtru[2][2] * matrice[j + (latime + 2) + 1] )) / Factor ; 

                                // In caz de e mai mic ca 0 rezultatul il voi face 0 .
                                if(temp[j] < 0){
                                    temp[j] = 0;
                                }
                                // In caz de rezultatul e mai mare ca maximumGray , ii voi da valoarea maxima.
                                if(temp[j] > maximumGray){
                                    temp[j] = maximumGray;
                                }
                            }
                            statistica[rank] += inaltime; 
                            MPI_Send(temp + (latime+2), (latime+2)*(inaltime), MPI_INT, parent, 5, MPI_COMM_WORLD);
                        }else {
                            MPI_Send(matrice + (latime+2), (latime+2)*(inaltime), MPI_INT, parent, 5, MPI_COMM_WORLD);
                        }
                            

                        ////////////////////// Am terminat de modificat pixelii                       
                    }    
                



             }

        }

        

        

    }

    MPI_Finalize();
    return 0;
}
