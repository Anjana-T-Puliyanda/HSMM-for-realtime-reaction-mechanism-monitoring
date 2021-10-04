function Z=elnproduct(elnx,elny)
 Z=0;
  
        if isnan(elnx)||isnan(elny)
            Z=NaN;
        else
            Z=elnx+elny;
        end


end