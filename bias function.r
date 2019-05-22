#Bias function for lme and nlme model fits (not applicable to models fitted using lme4)

#Re-written to incorporate mean % error and mean absolute % error
#Dave Auty - 04/12/2012

bias=function(pred,obs) 
  	{
		if (length(pred)!=length(obs)) warning ("length has to be equal")

	bias=array(dim=5)
		
 bias[1]=mean(obs-pred)                   #mean error
	bias[2]=mean(abs(obs-pred))              #mean absolute error (model accuracy)
	bias[3]=sqrt(mean((obs-pred)^2))         #RMSE (average magnitude of errors)
 bias[4]=mean((obs-pred)/pred)*100        #mean percent error (new - different to [5])
	bias[5]=mean(abs(obs-pred)/pred)*100     #mean absolute percent error (Eq.[9] in Parresol 1999, but confusingly named as 'mean percent standard error (S%))
    		
 names(bias)=c("E", "|E|", "RMSE", "E%", "|E|%")
		 
		return (bias)
			}
