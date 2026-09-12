cd /p/scratch/chwb03-2/hwb03f/calculations/monolayer_graphene/lattice12x12/V9.1_U3.0


perl -pe 's/^/0  /' ./processing_1el/GF_single_el_K1_0_K2_0_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func0.txt 
perl -pe 's/^/1  /' ./processing_1el/GF_single_el_K1_0_K2_1_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func1.txt 
perl -pe 's/^/2  /' ./processing_1el/GF_single_el_K1_0_K2_2_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func2.txt 
perl -pe 's/^/3  /' ./processing_1el/GF_single_el_K1_0_K2_3_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func3.txt 
perl -pe 's/^/4  /' ./processing_1el/GF_single_el_K1_0_K2_4_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func4.txt 
perl -pe 's/^/5  /' ./processing_1el/GF_single_el_K1_0_K2_5_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func5.txt 
perl -pe 's/^/6  /' ./processing_1el/GF_single_el_K1_0_K2_6_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func6.txt 
perl -pe 's/^/7  /' ./processing_1el/GF_single_el_K1_2_K2_7_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func7.txt 
perl -pe 's/^/8  /' ./processing_1el/GF_single_el_K1_4_K2_8_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func8.txt 
perl -pe 's/^/9  /' ./processing_1el/GF_single_el_K1_3_K2_6_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func9.txt 
perl -pe 's/^/10  /' ./processing_1el/GF_single_el_K1_2_K2_4_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func10.txt 
perl -pe 's/^/11  /' ./processing_1el/GF_single_el_K1_1_K2_2_S1_0_S2_0_L1_0_L2_0/logs_an_cont/averaged_analytical_cont_spectral_func_final.txt  > ./sp_func11.txt 




perl -ne  'print if $.> 2' ./sp_func0.txt   >  ./sp_func0_1.txt 
perl -ne  'print if $.> 2' ./sp_func1.txt   >  ./sp_func1_1.txt 
perl -ne  'print if $.> 2' ./sp_func2.txt  >  ./sp_func2_1.txt 
perl -ne  'print if $.> 2' ./sp_func3.txt >  ./sp_func3_1.txt 
perl -ne  'print if $.> 2' ./sp_func4.txt >  ./sp_func4_1.txt 
perl -ne  'print if $.> 2' ./sp_func5.txt >  ./sp_func5_1.txt 
perl -ne  'print if $.> 2' ./sp_func6.txt >  ./sp_func6_1.txt 
perl -ne  'print if $.> 2' ./sp_func7.txt >  ./sp_func7_1.txt 
perl -ne  'print if $.> 2' ./sp_func8.txt >  ./sp_func8_1.txt 
perl -ne  'print if $.> 2' ./sp_func9.txt >  ./sp_func9_1.txt 
perl -ne  'print if $.> 2' ./sp_func10.txt >  ./sp_func10_1.txt 
perl -ne  'print if $.> 2' ./sp_func11.txt >  ./sp_func11_1.txt 


perl -ne  'print; END {print "\n"}' ./sp_func0_1.txt   >  ./sp_func0_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func1_1.txt   >  ./sp_func1_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func2_1.txt  >  ./sp_func2_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func3_1.txt >  ./sp_func3_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func4_1.txt >  ./sp_func4_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func5_1.txt >  ./sp_func5_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func6_1.txt >  ./sp_func6_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func7_1.txt >  ./sp_func7_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func8_1.txt >  ./sp_func8_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func9_1.txt >  ./sp_func9_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func10_1.txt >  ./sp_func10_2.txt 
perl -ne  'print; END {print "\n"}' ./sp_func11_1.txt >  ./sp_func11_2.txt 


cat ./sp_func0_2.txt ./sp_func1_2.txt   ./sp_func2_2.txt  ./sp_func3_2.txt  ./sp_func4_2.txt  ./sp_func5_2.txt  ./sp_func6_2.txt  ./sp_func7_2.txt  ./sp_func8_2.txt  ./sp_func9_2.txt  ./sp_func10_2.txt  ./sp_func11_2.txt   > all_func_1el.txt 
rm -rf ./sp_func*.txt
