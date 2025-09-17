using bestNet
using Test


@testset "betnet.jl" begin
    #Load the temporal and Genomic Data
    @time transNetworkInference(tempfile="./bestNet/data/time_real_data.csv",SNPfile="./betnet2.0/data/SNP_real_data.csv",Contactfile="",genomeSize=4411532, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_69.csv")
    
    # simulation data with/without contact network
    @time transNetworkInference(tempfile="./bestNet/data/TemporalData1_100.csv",SNPfile="./betnet2.0/data/SNP1_100.csv",Contactfile="",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_nonet.csv")
    @time transNetworkInference(tempfile="./bestNet/data/TemporalData1_100.csv",SNPfile="./betnet2.0/data/SNP1_100.csv",Contactfile="./betnet2.0/data/ContactProb1_100.csv",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_net.csv")
 
end
