function DELTA =CalcDelta(V,opts,k)
    if nargin < 3 || isempty(opts)
        opts = "default";
    end
    
    switch lower(string(opts))
        case "default"
            DELTA = rand(1,V)/V; 
        case "sin"
            DELTA=zeros(1,V);
            for i=1:V
                DELTA(i) = sin(k)/V;
            end
        case "max"
            DELTA=zeros(1,V);
            DELTA(randi(V))=1;
    end
end