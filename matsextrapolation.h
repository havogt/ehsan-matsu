template<typename Function>
double extrapolateMats(const int n,const Function &F,const int Mats,const bool extraLogarithmic = true)
{

            double result=0.0;

            if(n>Mats-1)
            {

                        double Fx1=F(Mats-1);
                        double Fx0=F(Mats-2);
                        double extra=0.0;
                        if (extraLogarithmic)
                        {
                                    extra=pow(10,log10(Fx0)+log10(Fx1/Fx0)*(n-(Mats-2))); //log-linear extrapolation
    
                        }
                        else
                                    extra=Fx0+(Fx1-Fx0)*(n-(Mats-2)); //linear extrapolation
                        result=extra;

            }

            else
            {
                        result=F(n) ;

            }

            return result;
}
