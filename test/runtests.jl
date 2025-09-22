using bestNet

using Test


#@testset "betnet.jl" begin
    #Load the temporal and Genomic Data
    @time transNetworkInference(tempfile="./bestNet/data/sim_temporal.csv",SNPfile="./bestNet/data/sim_snp.csv",Contactfile="",genomeSize=100000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_sim_50.csv")
    
    # simulation data with/without contact network
    @time transNetworkInference(tempfile="./bestNet/data/TemporalData1_100.csv",SNPfile="./betnet2.0/data/SNP1_100.csv",Contactfile="",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_nonet.csv")
    @time transNetworkInference(tempfile="./bestNet/data/TemporalData1_100.csv",SNPfile="./betnet2.0/data/SNP1_100.csv",Contactfile="./betnet2.0/data/ContactProb1_100.csv",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_net.csv")
 
#end
