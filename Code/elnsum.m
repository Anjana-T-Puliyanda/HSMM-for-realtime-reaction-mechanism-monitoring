function Z=elnsum(elnx,elny)
Z=0;


        if isnan(elnx)||isnan(elny)
            if isnan(elnx)
                Z=elny;
            else
                Z=elnx;
            end
        else
            if elnx>elny
                Z=elnx+eln(1+exp(elny-elnx));
            else
                Z=elny+eln(1+exp(elnx-elny));
            end
        end              


end