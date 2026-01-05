#include <Rcpp.h>
#include <cmath>       /* tgamma */
#include <limits>
#include <vector>       /* tgamma */
#include <unordered_set>
#include <string>
#include <unordered_map>
#include <set>
#include <utility>
//#include "model.hpp"

using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
// http://www.rcpp.org/
// http://adv-r.had.co.nz/Rcpp.html
// http://gallery.rcpp.org/


//https://www.boost.org/doc/libs/1_58_0/libs/math/doc/html/math_toolkit/dist_ref/dists/students_t_dist.html
//https://www.boost.org/doc/libs/1_58_0/libs/math/doc/html/math_toolkit/sf_beta/beta_function.html
//double t_bf_sd= 0.12;//0.08397286;//这个两个参数基本能够让两个分布算出来的llv值接近

//double t_ratio_sd=0.16;//0.3181398; //调整t_ratio的理论计算方差，更多的根据ratio来进行划分

//vector<double> p_cnts{0.1, 0.15, 0.25, 0.15, 0.1, 0.1, 0.1 ,0.05};//pn_t先验概率，这组参数会导致整体的ploidy偏低
//vector<double> p_cnts{0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1 ,0.1};//pn_t先验概率
//vector<double> p_cnts{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ,0.1};//pn_t先验概率

/*
double logit(double x){
  double y=1.6*x+0.1;
  return log( y/(1-y) );
}



// Normal PDF(x) = exp(-x*x/2)/{sigma * sqrt(2 * Pi) }

double dnorm(const double & x){
  double Pi= 3.14159265358979323846 ;

   //double a=(1.0/(double)pow(2 * Pi, 0.5)) * exp(-0.5 * x * x);

  double b=0-log( (double)pow(2*Pi,0.5) )+(-0.5*x*x);

  //return log(a);
  return b;
}

NumericVector dnormr(NumericVector x) {
  NumericVector y;
  //cout<<dnorm( x.at(0) )<<endl;

  return y;
}

void test_log(NumericVector x) {
  //init_log_table();
  //std::cout<<fast_log2(100)<<" "<<log2(100);

}

*/

//https://www.appsloveworld.com/cplus/100/39/very-fast-approximate-logarithm-natural-log-function-in-c
//constexpr int LogPrecisionLevel = 20;
//constexpr int LogTableSize = 1 << LogPrecisionLevel;
//double* log_table=new double[LogTableSize];// [[Rcpp::export]]

constexpr int LogPrecisionLevel = 20;
constexpr int LogTableSize = 1 << LogPrecisionLevel;
static double* log_table=new double[LogTableSize];




//[[Rcpp::export(rng = false)]]
void init_log_table() {
  for (int i = 0; i < LogTableSize; i++) {
    log_table[i] = log2(1 + (double)i / LogTableSize);
  }
  //cout<<LogTableSize<<endl;

}



inline double fast_log2(double x) { // x>0
  long long t = *(long long*)&x;
  int exp = (t >> 52) - 0x3ff;
  int mantissa = (t >> (52 - LogPrecisionLevel)) & (LogTableSize - 1);
  return exp + log_table[mantissa];
}





/*** R
init_log_table();
*/

double change_log_v=log2(exp(1));

double v=5;

long double normlize_constant=(sqrt(v) *tgamma(v/2)*tgamma(0.5)/tgamma(v/2+ 0.5) );

double inline dt(double x){
  //long double dt_v=pow( (v / (v + x*x)),(1+v)/2) / (sqrt(v) *tgamma(v/2)*tgamma(0.5)/tgamma(v/2+ 0.5) );

  long double dt_v=pow( (v / (v + x*x)),(1+v)/2) / normlize_constant;

  return log2(dt_v);
}


double normlize_constant_log=log2(sqrt(v) *tgamma(v/2)*tgamma(0.5)/tgamma(v/2+ 0.5) );

inline double dt2(double x){
  //long double dt_v=pow( (v / (v + x*x)),(1+v)/2) / (sqrt(v) *tgamma(v/2)*tgamma(0.5)/tgamma(v/2+ 0.5) ) ;

  //long double dt_v=pow( (v / (v + x*x)),(1+v)/2) / normlize_constant ;
  double dt_v=((1+v)/2)*fast_log2( (v / (v + x*x)) )-  normlize_constant_log ;

  return dt_v;
}

double v2=5;

double normlize_constant_log2=log2(sqrt(v2) *tgamma(v2/2)*tgamma(0.5)/tgamma(v2/2+ 0.5) );

inline double dt3(double x){
  //long double dt_v=pow( (v / (v + x*x)),(1+v)/2) / (sqrt(v) *tgamma(v/2)*tgamma(0.5)/tgamma(v/2+ 0.5) ) ;

  //long double dt_v=pow( (v / (v + x*x)),(1+v)/2) / normlize_constant ;
  double dt_v=((1+v2)/2)*fast_log2( (v2 / (v2 + x*x)) )-  normlize_constant_log2 ;

  return dt_v;
}

//fast_log2

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

  /*** R
# timesTwo(42)
*/



inline double baf_bayes_one_fit_cpp_v(vector<double> &ratio, vector<double> &bf,
                                    vector<int> &segs,vector<string> &zygosity_v,
                                    int start_i, int end_i,
                                    vector<int> &n_v, vector<int> &m_v,
                                    vector<double> &trunca_ratio_t_v,
                                    vector<double> &trunca_baf_t_v,

                                    vector<double> &ratio_v,
                                    vector<double> &baf_v,
                                    double ratio_sd_v,
                                    double bf_sd_v,
                                    double bf_sd_logit_v){

  int const_larget_pen=-1*100000;

  if(end_i-start_i+1<3){
    return(const_larget_pen);
  }

  double seg_len=(segs.at(end_i)-segs.at(start_i)+1)/1e6;

  //cout<<"seg_len="<<seg_len<<endl;
  if(seg_len<1){
    return(const_larget_pen);
  }


  int CNn=2;

  double best_llv=-1 * std::numeric_limits<double>::max();
  int best_CNt=0;
  int best_B=0;
  int best_A=0;

  int lastn=-1;
  double last_llv=-1 * std::numeric_limits<double>::max();


  for(int k=0; k<n_v.size(); ){

    int CNt=n_v.at(k);
    int B=m_v.at(k);


    //theorectical ratio
    //double exp_ratio=(( (B * cellularity) + ( 1 - cellularity) )+ ( ((CNt-B) * cellularity) + 1 * ( 1 - cellularity) ) )/2;

    //theorectical bf
    //double exp_baf = ( (B * cellularity) + ( 1 - cellularity) ) /
    //  ( (CNt * cellularity) + CNn * ( 1 - cellularity) +0.000001);

    double exp_ratio=ratio_v.at(k);
    double exp_baf=baf_v.at(k);

    double trunca_ratio_t_v_one=trunca_ratio_t_v.at(k);
    double trunca_baf_t_v_one=trunca_baf_t_v.at(k);

    double llv_one=0;
    for(int i=start_i;i<=end_i;i++){

      llv_one+=dt2( ((ratio.at(i) )-( exp_ratio) )/ratio_sd_v ) - trunca_ratio_t_v_one;


      //llv_one+=dt( ( (ratio.at(i) )- ( exp_ratio) )/ratio_sd_v ) - trunca_norma_ratio_v_one;
      //cout<<"llv_one ratio: "<<llv_one<<endl;

      //if(bf.at(i)!=0){
      if( ( zygosity_v.at(i)=="het") && (CNt>0)  ){
        llv_one+=dt3( (bf.at(i)-exp_baf)/bf_sd_v ) - trunca_baf_t_v_one;//truncated normal dis

        //cout<<"llv_one bf: "<<llv_one<<endl;
        //llv_one+=dnorm( (logit( 1.9999998*bf.at(i)+0.00000001 ) - logit( 1.9999998*exp_baf+0.00000001) )/bf_sd_logit_v ); //baf*2在0-1。logit-normal dis
        //llv_one+=dnorm( (logit( bf.at(i) ) - logit( exp_baf ) )/bf_sd_logit_v )-trunca_baf_t_v_one; //baf*2在0-1。logit-normal dis

      }

    }

    //llv_one=llv_one/log2(exp(1));
    llv_one=llv_one/change_log_v;

    //影响了整体ploidy的计算？
    //llv_one=llv_one+log(p_cnts.at(CNt) );//pn_t先验概率
    //
    //std::cout<<llv_one<<std::endl;
    //std::cout<<best_llv<<std::endl;

    if(llv_one>best_llv){

      best_CNt=CNt;
      best_B=B;
      best_A=CNt-B;

      best_llv=llv_one;

    }

    if( ( CNt==lastn)&&(llv_one<last_llv) ){


      //while( (++k<= (n_v.size()-1)) && ( n_v.at(k)==lastn) ){
        //std::cout<<"k="<<k<<std::endl;


      //};
      //continue;
    }


    lastn=CNt;
    last_llv=llv_one;

    k++;

  }

  //NumericVector re;

  //re.push_back(best_CNt);
  //re.push_back(best_A);
  //re.push_back(best_B);
  //re.push_back(best_llv);

  return(best_llv);
}



//[[Rcpp::export(rng = false)]]
List baf_bayes_fit_cpp(NumericVector ratio, NumericVector bf, NumericVector seg,
                       StringVector zygosity,
                                    NumericVector n_exp, NumericVector m_exp,
                                    NumericVector baf_exp, NumericVector ratio_exp,
                                    NumericVector trunca_ratio_t, NumericVector trunca_baf_t,

                                    NumericVector potential_pos,
                                    NumericVector beta,
                                    NumericVector ratio_sd,
                                    NumericVector bf_sd,
                                    NumericVector t_bf_sd_logit
                                    ){
  //init_log_table();


  double ratio_sd_v=ratio_sd.at(0);

  double bf_sd_v=bf_sd.at(0);

  double t_bf_sd_logit_v=t_bf_sd_logit.at(0);


  vector<int> n_exp_v;
  for(int i=0;i<n_exp.size();i++){n_exp_v.push_back(n_exp.at(i));}

  vector<int> m_exp_v;
  for(int i=0;i<m_exp.size();i++){m_exp_v.push_back(m_exp.at(i));}

  vector<double> baf_exp_v;
  for(int i=0;i<baf_exp.size();i++){baf_exp_v.push_back(baf_exp.at(i));}

  vector<double> ratio_exp_v;
  for(int i=0;i<ratio_exp.size();i++){ratio_exp_v.push_back(ratio_exp.at(i));}

  //NumericVector trunca_norma_ratio, NumericVector trunca_baf_ratio,
  //vector<double> baf_exp_v;
  //for(int i=0;i<baf_exp.size();i++){baf_exp_v.push_back(baf_exp.at(i));}

  vector<double> trunca_ratio_t_v;
  for(int i=0;i<trunca_ratio_t.size();i++){trunca_ratio_t_v.push_back(trunca_ratio_t.at(i));}

  vector<double> trunca_baf_t_v;
  for(int i=0;i<trunca_baf_t.size();i++){trunca_baf_t_v.push_back(trunca_baf_t.at(i));}

  //vector<double> trunca_baf_logistic_v;
  //for(int i=0;i<trunca_baf_logistic.size();i++){trunca_baf_logistic_v.push_back(trunca_baf_logistic.at(i));}

  double beta_v=beta.at(0);

  vector<double> ratio_v;
  for(int i=0;i<ratio.size();i++){ratio_v.push_back(ratio.at(i));}

  vector<double> bf_v;
  for(int i=0;i<bf.size();i++){bf_v.push_back(bf.at(i));}

  vector<int> seg_v;
  for(int i=0;i<seg.size();i++){seg_v.push_back(seg.at(i));}

  vector<std::string> zygosity_v;
  for(int i=0;i<zygosity.size();i++){zygosity_v.push_back(Rcpp::as<std::string>(zygosity.at(i)) );}



  unordered_set<int> potential_pos_v;
  for(int i=0;i<potential_pos.size();i++){potential_pos_v.insert(potential_pos.at(i));}


  if(potential_pos_v.find(1) != potential_pos_v.end()){

    //cout<<"x"<<endl;
  }


  if(potential_pos_v.find(1) != potential_pos_v.end()){

    //cout<<"y"<<endl;
  }


  // for(int i=0;i<=potential_pos_v.size();i++){
  //   for(int j=0;j<=potential_pos_v.size();j++){
  //
  //   }
  //
  // }


  //vector<double> fv;fv.reserve(ratio_v.size()+1 );fv.push_back(3.087308e-12);fv.push_back( 1.000000e+05);
  //vector<int> cp;fv.reserve(ratio_v.size()+1 );cp.push_back(1);cp.push_back(1);

  vector<double> fv;for( int k=0;k<ratio_v.size();k++){fv.push_back(-1); }
  vector<int> cp;for( int k=0;k<ratio_v.size();k++){cp.push_back(-1); }

  fv[0]=(3.087308e-12);fv[1]=( 1.000000e+05);
  cp[0]=(1);cp[1]=(1);
  //cout<<"loop"<<endl;


  set<int> potential_pos_order_v;
  for(auto xx:potential_pos_v){potential_pos_order_v.insert(xx);}

  set<int> potential_pos_order_vv;
  for(auto xx:potential_pos_v){potential_pos_order_vv.insert(xx);}

  for( int i: potential_pos_order_vv){

    if(i<3){continue;}

  //for(int i=3;i<=ratio_v.size();i++){


    //if(i % 500==0){
      //cout<<"i="<<i<<endl;
    //}
    //cout<<"i="<<i<<" ratio size= "<<ratio_v.size()<<endl;


    if (potential_pos_v.find(i) == potential_pos_v.end() ){
      //cp.push_back(cp[cp.size()-1]);
      //fv.push_back(fv[fv.size()-1]+0.1);
      //cp[i-1]=cp[i-2];//cout<<"x="<<endl;
      //fv[i-1]=fv[i-2]+0.1;//cout<<"y="<<endl;

      continue;
    }



    int index_min_j=10000000;
    double min_j=10000000;


    unordered_map<int, double> llv_j;
    //cout<<"start"<<endl;
    for( int j: potential_pos_order_v){

      if(j>i-1){break;}

      //if (potential_pos_v.find(j) == potential_pos_v.end()){
      //  continue;
      //}

    /*
    for(int j=1;j<=i-1;j++){


      if (potential_pos_v.find(j) == potential_pos_v.end()){
        continue;
      }
    */


     //cout<<"i="<<i<<"j="<<j<<endl;
      //cout<<"aa"<<endl;

      //cout<<"x "<<j+1-1<<" "<<i-1<<endl;


     //这里的乘以2，是因为这里所有的cost都统一x2，包括了beta

     double tmp_jj=(-2)*baf_bayes_one_fit_cpp_v( ratio_v,
                    bf_v, seg_v,
                    zygosity_v,
                    j+1-1, i-1,
                    n_exp_v, m_exp_v,
                    trunca_ratio_t_v,trunca_baf_t_v,
                    ratio_exp_v, baf_exp_v,
                    ratio_sd_v, bf_sd_v, t_bf_sd_logit_v);

    llv_j[(j+1-1)*100000+ i-1] = tmp_jj;


    double tmp_j=fv.at(j-1)+tmp_jj+beta_v;

      //cout<<"aaa"<<endl;

      // if(i==18){
      //   cout<<"tmp_j_1="<<j<<endl;
      //
      //   cout<<"tmp_j_1="<<seg_v.at(i-1)<<endl;
      //   cout<<"tmp_j_1="<<seg_v.at(j+1-1)<<endl;
      //
      //   cout<<"tmp_j_1="<<tmp_j<<endl;
      //   cout<<"tmp_j_1="<<fv.at(j-1)<<endl;
      //
      //   cout<<"tmp_j_2="<<(-2)*baf_bayes_one_fit_cpp_v( ratio_v, bf_v, seg_v,
      //              j+1-1, i-1,
      //              n_exp_v,m_exp_v,ratio_exp_v,baf_exp_v )<<endl;
      //
      // }


      //double tmp_j=1;


      //cout<<"tmp_j_3="<<beta_v<<endl;


      if(tmp_j<min_j){
        index_min_j=j;
        min_j=tmp_j;

      }

    }

    if(index_min_j==10000000){
      //cp.push_back(cp[cp.size()-1]);
      //fv.push_back(fv[fv.size()-1]+0.1);
      cp[i-1]=cp[i-2];//cout<<"xx="<<endl;
      fv[i-1]=fv[i-2]+0.1;//cout<<"yy="<<endl;

      //continue;

    }else{
      //cp.push_back(index_min_j);
      //fv.push_back(min_j);

      cp[i-1]=index_min_j;//cout<<"xxx="<<endl;
      fv[i-1]=min_j;//cout<<"yyy="<<endl;
      //continue;
    }


    //cout<<"pelt"<<endl;

    //PELT方法，数据量大的时候能够有一个比较明显的加速。

    vector<int> potential_pos_v_remove;

    unordered_set<int>::iterator potential_pos_iter = potential_pos_v.begin();

    potential_pos_iter++;

    if(ratio_v.size()>700){

      while (potential_pos_iter != potential_pos_v.end()){
        if( (*potential_pos_iter)>=i ){
          potential_pos_iter++;
          continue;
        }

        //cout<<"y "<<*potential_pos_iter+1-1<<" "<<i-1<<endl;
       /*
        if(fv[*potential_pos_iter-1]+
           (-2)*baf_bayes_one_fit_cpp_v( ratio_v, bf_v, seg_v,
            zygosity_v,*potential_pos_iter+1-1, i-1,n_exp_v, m_exp_v,
            trunca_ratio_t_v,trunca_baf_t_v,
            ratio_exp_v, baf_exp_v,
            ratio_sd_v, bf_sd_v, t_bf_sd_logit_v)-beta_v > fv[i-1] ){
        */
       if(fv[*potential_pos_iter-1]+
          llv_j[(*potential_pos_iter+1-1)*100000+ i-1]-beta_v > fv[i-1] ){


          //std::cout<<"tau: "<<*potential_pos_iter<<endl;


          //std::cout<<"i: "<<i<<endl;

          //std::cout<<"f(tau): "<<fv[*potential_pos_iter-1]<<std::endl;
          //std::cout<<"cost: "<<(-2)*baf_bayes_one_fit_cpp_v( ratio_v, bf_v, seg_v,
          //                zygosity_v,*potential_pos_iter+1-1, i-1,n_exp_v, m_exp_v,
          //                ratio_exp_v, baf_exp_v )<<std::endl;

          //std::cout<<"beta: "<<beta_v<<endl;

          //std::cout<<"f(t): "<<fv[i-1]<<endl;

          //std::cout<<"remove: "<< (*potential_pos_iter)<<std::endl;

		      potential_pos_v_remove.push_back(*potential_pos_iter);
          //potential_pos_v.erase(*potential_pos_iter);

        }

        potential_pos_iter++;

      }

	    for(int k=0;k<potential_pos_v_remove.size();k++){
	      potential_pos_v.erase(potential_pos_v_remove.at(k));
	      potential_pos_order_v.erase(potential_pos_v_remove.at(k));
	    }

     }

    //cout<<"end"<<endl;


    // if(index_min_j==1 && i==23){
    //   cout<<fv.at(index_min_j-1)<<endl;
    //   cout<<(-2)*baf_bayes_one_fit_cpp_v( ratio_v, bf_v, seg_v,
    //          index_min_j+1-1, i-1,
    //          n_exp,m_exp,ratio_exp,baf_exp )<<endl;
    //
    //   cout<<beta_v<<endl;
    //
    //   cout<<index_min_j<<endl;
    //
    //   cout<<min_j<<endl;
    //
    // }


  }

  //cout<<"output"<<endl;


  NumericVector fv_r(fv.size());
  NumericVector cp_r(cp.size());

  for(int i=0;i<fv.size();i++){fv_r[i]=(fv.at(i));}
  for(int i=0;i<cp.size();i++){cp_r[i]=(cp.at(i));}

  List best_seg = List::create(Named("fv") = fv_r , _["cp"] = cp_r);

  return best_seg;

}




// [[Rcpp::export]]
NumericVector baf_bayes_one_fit_cpp(NumericVector ratio, NumericVector bf, StringVector zygosity,
                                    NumericVector n_v, NumericVector m_v,
                                    NumericVector ratio_v,NumericVector baf_v,
                                    NumericVector trunca_ratio_t_v,NumericVector trunca_baf_t_v,
                                    NumericVector ratio_sd,
                                    NumericVector bf_sd,
                                    NumericVector bf_sd_logit){

  //init_log_table();


  int CNn=2;

  //double cellularity=t_cellularity.at(0);

  double best_llv=-1 * std::numeric_limits<double>::max();
  int best_CNt=0;
  int best_B=0;
  int best_A=0;

  for(int k=0; k<n_v.size(); k++){

    int CNt=n_v.at(k);
    int B=m_v.at(k);

    //
    //cout<<" "<<n_v.at(k)<<" "<<m_v.at(k)<<" "<<ratio_v.at(k)<<" "<<baf_v.at(k)<<" " <<endl;
    //theorectical ratio;
    //double exp_ratio=(( (B * cellularity) + ( 1 - cellularity) )+ ( ((CNt-B) * cellularity) + 1 * ( 1 - cellularity) ) )/2;
    //theorectical bf
    //double exp_baf = ( (B * cellularity) + ( 1 - cellularity) ) /
    //  ( (CNt * cellularity) + CNn * ( 1 - cellularity) +0.000001);


    double exp_ratio=ratio_v.at(k);

    double exp_baf=baf_v.at(k);

    double trunca_ratio_t_v_one=trunca_ratio_t_v.at(k);

    double trunca_baf_t_v_one=trunca_baf_t_v.at(k);


    //cout<<" trunca_baf_t_v_one: "<<trunca_baf_t_v_one<<" ";

    //double llv_one_sum=0;
    double llv_one=0;

    for(int i=0;i<ratio.size();i++){
      //double llv_one;
      llv_one+=dt2( ( (ratio.at(i) )-( exp_ratio) )/ratio_sd.at(0) )-trunca_ratio_t_v_one;
      //double llv_one=dnorm( ((ratio.at(i) )-( exp_ratio) )/ratio_sd.at(0) )-trunca_norma_ratio_v_one;

      //llv_one+=dt( (log(ratio.at(i) )-log( exp_ratio) )/ratio_sd.at(0) )-trunca_norma_ratio_v_one;
      //cout<<"llv_one ratio: "<<llv_one<<" "<<(ratio.at(i))<<" "<<( exp_ratio)<<" "<<i <<" "<<ratio.size()<<endl;
      //cout<<"llv_one ratio: "<<llv_one<<" "<<((ratio.at(i) )-( exp_ratio) )/ratio_sd.at(0)<<" "<<trunca_ratio_t_v_one<<" "<<i <<" "<<ratio.size()<<" "<<k<<endl;

      //llv_one_sum=llv_one_sum+llv_one;

      //if(bf.at(i)!=0){
      //if(zygosity.at(i)=="het"){
      //CNt=0的话求baf没意义啊
      if( ( zygosity.at(i)=="het") && (CNt>0)  ){

        llv_one+=dt2( (bf.at(i)-exp_baf)/bf_sd.at(0) )-trunca_baf_t_v_one;// truncated t dis
        //llv_one=dnorm( (logit( 1.997*bf.at(i)+0.001 )- logit( 1.997*exp_baf+0.001) )/bf_sd_logit.at(0) );//baf*2在0-1。logit-normal
        //llv_one+=dnorm( (logit( bf.at(i) ) - logit( exp_baf ) )/bf_sd_logit.at(0) )-trunca_baf_t_v_one; //baf*2在0-1。logit-normal dis

        //cout<<"llv_one bf: "<<llv_one<<" "<<(bf.at(i))<<" "<<( exp_baf)<<" "<<bf_sd.at(0)<<" "<<i<<endl;//break;
        //cout<<"llv_one bf: "<<llv_one<<" "<<(logit( bf.at(i) )- logit( exp_baf ) )/bf_sd_logit.at(0) << " "<<trunca_baf_t_v_one<<" "<<i<<endl;break;


        //llv_one_sum=llv_one_sum+llv_one;

        //1.998而不是2是为了+0.001，这样bf为0时不会导致log(0)，bf=0.5时，不会导致log负值
      }
    }

    llv_one=llv_one/log2(exp(1));


    //llv_one=llv_one+log(p_cnts.at(CNt) );//pn_t先验概率
    //std::cout<<llv_one<<std::endl;

    if(llv_one>best_llv){

      best_CNt=CNt;
      best_B=B;
      best_A=CNt-B;

      best_llv=llv_one;

    }


  }

  NumericVector re;

  re.push_back(best_CNt);
  re.push_back(best_A);
  re.push_back(best_B);
  re.push_back(best_llv);

  return(re);
}




