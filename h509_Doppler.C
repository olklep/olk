#define h509_Doppler_cxx
#include "h509_Doppler.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "TH1F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"

	void h509_Doppler::Loop(Long64_t numEntries)
	{	
  	//   In a ROOT session, you can do:
	//      Root > .L h509_Doppler.C
  	//      Root > h509_Doppler t
  	//      Root > t.GetEntry(12); // Fill t data members with entry number 12
  	//      Root > t.Show();       // Show values of entry 12
  	//      Root > t.Show(16);     // Read and show values of entry 16
  	//      Root > t.Loop();       // Loop on all entries
  	//

  	//    This is the loop skeleton where:
  	//    jentry is the global entry number in the chain
  	//    ientry is the entry number in the current Tree
  	//    Note that the argument to GetEntry must be:
  	//    jentry for TChain::GetEntry
  	//    ientry for TTree::GetEntry and TBranch::GetEntry
  	//
  	//    To read only selected branches, Insert statements like:
  	// 	METHOD1:
  	//    fChain->SetBranchStatus("*",0);  // disable all branches
  	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
  	// METHOD2: replace line
  	//    fChain->GetEntry(jentry);       //read all branches
  	//	by  b_branchname->GetEntry(ientry); //read only this branch
  		

	if (fChain == 0) return;

	 	 // Define Histos
	 	 //This plots all 162 crystals of the crystal ball without any conditions 		applied 
 	 	Int_t i;
 	 	TH1F *xbe[162];
 	 	for (i = 0; i < 162; i++) 
	{
 	 	  char name[512];
 	 	  sprintf(name, "xbe_%03d", i+1);
 	 	  xbe[i] = new TH1F(name, "E_SYNC_XB", 8192, 0, 8192);
  	 	
 	 		}	
  	
  	TH1F *Energy_Sum = new TH1F("Energy_sum", "XB_Sume", 100000, 0, 													100000);//physical sum
	TH1F *Energy_sum_histo = new TH1F("sum_histo", "sum_histo", 8192, 											0.,8192.);
	TH2F *Energy_sum_multi = new TH2F("sum_multi", "sum_multi",100,0.,100., 									200000, 0.,20000.);
	//TH2F *h_PID = new TH2F("PID","PID",700,0.,35.,200,0.,5.); 
	  //TH2F *h_PID_cut = new TH2F("PID_cut","PID_cut",700,0.,35.,200,0.,5.);
	  //TH2F *h_Histo = new TH2F("Something","something"100,0.,180, 						100,0.,10);
	  //plots beta from the land02
	  TH1F* h_beta = new TH1F("beta","beta",1000,0.,1.);
	  TH1F* h_fra = new TH1F("fra","fra",100,25.,35.);
	
	  TH1F* h_nb_cluster = new TH1F("nb_cluster","nb_cluster",162,0.,162.);
	  TH1F* h_cluster_size = new 										 				  TH1F("cluster_size","cluster_size",162,0.,162.);
	  TH2F* h_nb_cluster_multi = new 	              TH2F("nb_cluster_multiplicity","nb_cluster_multiplicity",162,0.,162.,162,0.,162.	);
  
	  //
	  TH1F *xbe_cluster[162];
	  for (i = 0; i < 162; i++) 
	    {
	      char name[512];
	      sprintf(name, "Energy_cluster_%03d", i+1);
	      xbe_cluster[i] = new TH1F(name, "E_cluster", 8192, 0., 8192.);
	      
	    }


    //***********************************************************
    //making Doppler correction 
	
	  ifstream infile;
	  infile.open("Doppler_Correction.txt");

	  double cosine_theta1[162];

	  if(infile.is_open())
	    {
	      while(!infile.eof())
	        {
		  int temp_nb;
	  
		  infile>>temp_nb;
		  infile>>cosine_theta1[temp_nb];
		}
	      cout<<"Doppler parameters loaded"<<endl;
	    }
	  else cout<<"E> input file can't be opened ! "<<endl;

	  infile.close();
  


	  //making acceptance of the gamma rays
	  //*********************************************************

	  //ifstream infile;
	  infile.open("Gamma_Acceptance.txt");

	  double cosine_theta2[162],theta[162];

	  if(infile.is_open())
	    {
	      while(!infile.eof())
	        {
		  int temp_nb;
	  
		  infile>>temp_nb;
		  infile>>cosine_theta2[temp_nb];
		  infile>>theta[temp_nb];
		}
	      cout<<"Gamma_Acceptance loaded"<<endl;
	    }
	  else cout<<"E> input file can't be opened ! "<<endl;
	
	  infile.close();



	  
	  //*************************************************************
  //reading neighbors crystals
  //ifstream infile;
	  infile.open("crystal.txt");
	
	  int neighbour[162][6];

	  if(infile.is_open())
	    {
	      while(!infile.eof())
	        {
		  int temp_nb;
		  double temp_dump;
		  
		  infile >> temp_nb;
		  infile >> temp_dump;
		  infile >> temp_dump;
		  infile >> temp_dump;
		  for(int i=0;i<6;i++)
		    {
		      int temp;
		      infile >> temp;
			  neighbour[temp_nb-1][i]=temp-1;
		    	}
		}
    	  cout<<"Neighbour loaded"<<endl;
    	}
  	else cout<<"E> input file can't be opened ! "<<endl;

  	infile.close();
  	/*
  	  for(int i=0;i<162;i++)
  	  {
  	    cout<<"crystal n"<<i<<" neighbour :";
  	    for(int j=0;j<6;j++)
		cout<<neighbour[i][j]<<" / ";
  	    cout<<endl;
  	  }
  

*/
  //***********************************************
  	//Event Loop
  	Long64_t nentries =0;
  	if(numEntries==0)//this we add when in loop we want to see not all events 	//but just a small statistics
	    nentries = fChain->GetEntriesFast();
	  else
	    nentries=numEntries;

	
	  cout << nentries << endl;
	
	 	 Int_t nbytes = 0, nb = 0;
	 	 //  Long64_t nbytes = 0, nb = 0
	 	 for (Long64_t jentry=0; jentry<nentries;jentry++) {
	 	   Long64_t ientry = LoadTree(jentry);
		
	 	   //   cout << jentry << endl;
	
			    if (ientry < 0) break;
	 	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	 	    //added as atrigger cut condition, that is good beam + NTF(Tpat&2)
	 	   //as trigger condition on good beam +CB sum,	when are summed only 			//those crystals which catch a signals in the crystal ball	
	    	//if(Tpat&8)
    


		//if(Tpat&16)

	    // {//--


	   
	   	  //this is PID we have to have a look in land02,identification particles
		     	//h_PID->Fill(Inz,Inaoverz);
   		 //taken from land02
   		 h_beta->Fill(Inbeta);
   		     h_fra->Fill(fra_A);  
  	 	 //take cut conditions only on Ar32
    



			//--
			if//((p_mul>1) 
	  	    (fra_A<32.5 && fra_A>31.5)
	  	    {

	
      


			//h_PID_cut->Fill(Inz,Inaoverz);




      
			// Filling histo 
			double Sume = 0;
			double Xbmule = 0; //multiplicity for all crystal ball 
			bool valid_event = true; 
			//cout << "Xbmul  " << Xbmul << endl;
			//for (i = 0; i < (int)Xbmul; i++){
			double Energy_crystal[162];
			for(int k =0;k<162;k++)
		  




		  Energy_crystal[k]=0.;
		bool crystal_todo[163];//163--
		for(int k=0;k<163;k++)//163--
		  crystal_todo[k]=false;

		for( i=0 ; i<Xbmul;i++)
		  {
		  
		    int det_number = Xbi[i]-1;
		    double Energy = Xbe[i];
		    //cout << "Xbe" << Xbe << endl;
		    if (Energy>0.) 
	      {
				Xbmule++;//real multiplicity 
		      }
	 
		   
	 //double gamma_a[163] = 1-sqrt(1-(inbeta^2));
		    //calculation for the gamma acceptance to see in which area of the 			//crystall ball gamma is more 
	
		    double temp_gamma_a = TMath::Sqrt(1.-Inbeta*Inbeta);
		    double gamma_a = 1./temp_gamma_a;
		  
	  //dopler correction is done by using the beta from land02
		  	  double Energy2 = Energy*gamma_a*(1.-								Inbeta*cosine_theta1[det_number]);
		    
		    	    
			    xbe[det_number]->Fill(Energy2);
			   // Energy_crystal[det_number]=Energy
			   Energy_crystal[det_number]=Energy2;// bilo 
			    crystal_todo[det_number]=true;
				Sume += Energy2;
		
			    if(det_number>162 || det_number<0) valid_event=false;
				  }
		

		/*int neighbour_1 = neighbour[det_number][0];
		  int neighbour_2 = neighbour[det_number][1];
		  int neighbour_3 = neighbour[det_number][2];
		  int neighbour_4 = neighbour[det_number][3];
		  int neighbour_5 = neighbour[det_number][4];
		  int neighbour_6 = neighbour[det_number][5];
	
		    
		
		  int E_neighbour_1 = Energy_crystal[neighbour_1];
		  int E_neighbour_2 = Energy_crystal[neighbour_2];
		  int E_neighbour_3 = Energy_crystal[neighbour_3];
		  int E_neighbour_4 = Energy_crystal[neighbour_4];
		  int E_neighbour_5 = Energy_crystal[neighbour_5];
	  	  int E_neighbour_6 = Energy_crystal[neighbour_6];
		*/
		
	int det_number = Xbi[0]-1;
			if(det_number<0 && det_number>162)
			  {
			    cout<<"det_number wrong"<<endl;
			    exit(-1);
			  }



			//cout<<"----------------------------"<<endl<<"EVENT : 					"<<jentry<<endl<<"--------------------------------"<<endl<<"start 		point:"<<det_number<<endl ;
			bool finish=false;
			//clustering,looking for the high energy
			int tab_center_cluster[162];
			double tab_E_total_center_cluster[162];
			int tab_cluster_size[162]; //how many neighbors crystals close to 				//each other had a signal

		for(int k=0;k<162;k++)
		  {
		    tab_center_cluster[k]=-1;
		    tab_E_total_center_cluster[k]=0.;
		    tab_cluster_size[k]=0;

		    // cout<<"crystal n"<<k<<" = "<<crystal_todo[k]<<endl;
		  }

		int nb_of_cluster=0;

	

		while(finish!=true)
		  {
		    // cout<<"new cluster"<<endl;
		    int close_neighbour[6];
		    double E_close_neighbour[6];
		    crystal_todo[det_number]=false;

		    bool center_cluster_found=false;
		    int center_cluster_index;
		    double E_total_cluster=0.;
		    int cluster_size=0;

		    while(center_cluster_found!=true)
		      {
			int temp_center_cluster_index=det_number;
			double temp_E_center_cluster=Energy_crystal[det_number];
			double temp_total_E_cluster =0.;
		
			for (int j=0;j<6;j++)
			  {

			    close_neighbour[j]=neighbour[det_number][j];
			    //cout<<"close_neighbour["<<j<<"] = "<<close_neighbour[j]<<" / ";
	
			    if(close_neighbour[j] > 161)
			      {
					E_close_neighbour[j]=-1.;
			      }
			    else
			      {
				E_close_neighbour[j]=Energy_crystal[close_neighbour[j]];
			      }
			    crystal_todo[close_neighbour[j]]=false;
		      	
			  }
			//cout<<endl<<"setting done"<<endl;


			for(int k=0;k<6;k++)
			  {
			    if(E_close_neighbour[k]> temp_E_center_cluster)
			      {
				temp_E_center_cluster=E_close_neighbour[k];
				temp_center_cluster_index=close_neighbour[k];
			      }
			  }

			if(temp_center_cluster_index==det_number)
			  {
			    center_cluster_found=true;
			    center_cluster_index=temp_center_cluster_index;
			    // cout<<"center_cluster_index:"<<center_cluster_index;
			    E_total_cluster=Energy_crystal[center_cluster_index];
			    for(int k=0;k<6;k++)
			      {
				//cout<<" / "<<center_cluster_index;
				int temp_index = neighbour[center_cluster_index][k];
				if(temp_index>161)
				  continue;
				else
				  {
				    cluster_size++;
				    E_total_cluster+=Energy_crystal[temp_index];
				    crystal_todo[temp_index]=false;
				  }
			      }
			    // cout<<endl<<"cluster done"<<endl;
			    bool next_found=false;
			    int next_index=0;
			    while(next_found!=true)
			      {
				if(crystal_todo[next_index]==false)
				  next_index++;
				else
				  {
				    // cout<<"next one is:"<<next_index<<endl;
				    next_found=true;
				    det_number=next_index;
				  }
				if(next_index>=162)
				  {
				    finish=true;
				    next_found=true;
				  }  
			      }
			  }
			else
			  {
			    det_number=temp_center_cluster_index;
			    //cluster_size++; //--
			  }
		      }
		    
		    tab_cluster_size[nb_of_cluster]=cluster_size;
		    tab_E_total_center_cluster[nb_of_cluster]=E_total_cluster;
		    tab_center_cluster[nb_of_cluster]=center_cluster_index;
		    nb_of_cluster++;
		    //cout<<"nb of cluster:"<<nb_of_cluster<<" / "<<"current 			center:"<<center_cluster_index<<endl;
			 if(nb_of_cluster>162)
		      {
			cout<<"E> nb_of_cluster >162"<<endl;
			exit(-1);
		      }

		  }
	
		h_nb_cluster->Fill(nb_of_cluster);
		h_nb_cluster_multi->Fill(nb_of_cluster,Xbmule);
	
		for(int k=0;k<nb_of_cluster;k++)
		  {
		    	xbe_cluster[tab_center_cluster[k]]>Fill(tab_E_total_center_cluster[k]);
		    h_cluster_size->Fill(tab_cluster_size[k]);
		  }
	
	
 
		if(Xbmul>0 && valid_event)// && det_number>0 && det_number<163)
		  {
		    Energy_sum_multi->Fill(Xbmule,Sume);
		    Energy_Sum->Fill(Sume);//physical sum for all events 
		  }
	    
		//}//-------tpat
	  }

	   }
	  //added this from the Tpat conditions
  for(int i=0;i<162;i++)
    {
      //only aading histograms to see for example the calibration is working or 	not,to get resolution
	      Energy_sum_histo->Add(xbe[i]);
	    }
  
	  double Inbeta2 = h_beta->GetMean();//it is for the getting beta from 	land02 and save it
	  double temp_gamma_a2 = TMath::Sqrt(1.-Inbeta2*Inbeta2);
	  double gamma_a2 = 1./temp_gamma_a2;
	    
	  double b[162],df_lab_df_cm[162];
	  for(int i=0;i<162;i++)
	   	 {
	   	 //formula for calculation ratio laboratory system and central mass......
      	b[i] = (1.- Inbeta2*cosine_theta2[i])*(1.- Inbeta2*cosine_theta2[i]);
      	df_lab_df_cm[i] = gamma_a2*gamma_a2*b[i];
    	}
	//to see graph for the gamma acceptance calculated simply by using a 	formula.theta is from land02
	  TGraph* g_acceptance = new TGraph(162,theta,df_lab_df_cm);
	  //g_acceptance->Draw("AC*");

	  TFile* outfile = new Tfile(name_output.c_str(),"RECREATE");//to write down
 		root_out file
  
  	outfile->cd();
  	Energy_sum_multi->Write();
  	Energy_Sum->Write();//physical
  	Energy_sum_histo->Write();//just adding the histograms
  	//  h_PID->Write();
  	//h_PID_cut->Write(); 
  	g_acceptance->Write();
  	for (i=0;i<162;i++)
  	  {
  	    xbe[i]->Write();
  	    xbe_cluster[i]->Write();
  	  }
  	h_cluster_size->Write();
  	h_nb_cluster->Write();
  	h_nb_cluster_multi->Write();
  	h_fra->Write();
  	outfile→Close();//to save analyzed file with results.


}

