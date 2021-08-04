for rep=5:6
    for t=1:192
        if stim_seq_day1(rep,t)<0 %this item should get a small response
            if response(rep,t)==13
                response_accuracy(rep,t)=1;
            else
                response_accuracy(rep,t)=0;
            end
        else
            if response(rep,t)==14
                response_accuracy(rep,t)=1;
            else
                response_accuracy(rep,t)=0;
            end
        end
    end
end