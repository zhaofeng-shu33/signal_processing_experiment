function launcher()
global sequence_length;
sequence_length=500;
close all
global random_sequence;
random_sequence=2*(rand(1,sequence_length)>0.5)-1;
global step_size;
error_square_sequence=zeros([1,sequence_length]);
%reset signal channel
global system_input_sequence
lms_part();
rls_part();
    function serialize_omega_array(omega_array,row_offset)
        final_str=sprintf('%.4f',omega_array(row_offset,1));
        for i=2:11
            final_str=strcat(final_str,' & ',sprintf('%.4f',omega_array(row_offset,i)));
        end
        disp('convolution result; signal channel and filter')
        disp(conv([0.3,0.9,0.3],omega_array(row_offset,:)));
        fileID = fopen('coefficent_out.tex','w');
        fprintf(fileID,final_str);
    end
    function lms_part()
        step_size=0.05;
        sigma_seq=[0.01,0.05,0.1];
        color_seq=['r','b','g'];
        one_time_randn=randn([1,sequence_length+2]);
        figure
        for j=1:3
        sigma=sigma_seq(j);
        system_input_sequence=conv(random_sequence,...
        [0.3,0.9,0.3])+sigma*one_time_randn;
        %lms system initialization
        lms_system(0,'c');
        for i=11:sequence_length
           [~,error_square_sequence(i)]=lms_system(i);
        end       
        %plot one time error_square convergent curve
        plot(error_square_sequence,color_seq(j));
        hold on
        end
        legend(sprintf('\\sigma=%f',sigma_seq(1)),...
        sprintf('\\sigma=%f',sigma_seq(2)),...
        sprintf('\\sigma=%f',sigma_seq(3)))
        title('一次实验LMS误差平方收敛曲线')
        xlabel('数字时间')
        ylabel('误差平方')
        saveas(gcf,'lml_one_time.eps','epsc')



        figure
        sigma=sigma_seq(2);
        step_size_seq=[0.1,0.05,0.01];
        error_square_sequence_3=zeros([3,sequence_length]);
        omega_out_array=zeros([3,11]);
        given_repeat_times=20;
        for repeat_times=1:given_repeat_times
        random_sequence=2*(rand(1,sequence_length)>0.5)-1;
        system_input_sequence=conv(random_sequence,[0.3,0.9,0.3])+...
        sigma*randn([1,sequence_length+2]);    
            for j=1:3
                step_size=step_size_seq(j);
                %lms system initialization
                lms_system(0,'c');
                for i=11:sequence_length
                   [~,errot_tmp]=lms_system(i);
                   error_square_sequence_3(j,i)=error_square_sequence_3(j,i)+errot_tmp;
                end
                 %get filter parameters
                 [~,~,omega_out]=lms_system(0,1);
                 omega_out_array(j,:)=omega_out_array(j,:)+omega_out;                 
            end
        %plot one time error_square convergent curve
        end
        omega_out_array=omega_out_array/repeat_times;
        serialize_omega_array(omega_out_array,2);
        error_square_sequence_3=error_square_sequence_3/given_repeat_times;
        for j=1:3
        plot(error_square_sequence_3(j,:),color_seq(j));
        hold on
        end
        legend(sprintf('步长=%f',step_size_seq(1)),...
        sprintf('步长=%f',step_size_seq(2)),...
        sprintf('步长=%f',step_size_seq(3)))
        title('不同步长误差平方均值收敛曲线')
        xlabel('数字时间')
        ylabel('误差平方均值')
        saveas(gcf,'lml_many_times.eps','epsc')
    end

    function rls_part()
        global forget_parameter;
        forget_parameter = 1;
        sigma=0.05;
        step_size=0.01;
        error_square_sequence_2=zeros([2,sequence_length]);
        system_input_sequence_initial=conv(random_sequence,[0.3,0.9,0.3]);
%the main purpose is for its comparision with lms_system;
        for k=1:20
            one_time_randn=randn([1,sequence_length+2]);        
            system_input_sequence=system_input_sequence_initial+sigma*one_time_randn;
            rls_system(0,'c');
            for i=11:sequence_length
               [~,error_tmp]=rls_system(i);
               error_square_sequence_2(1,i)=error_square_sequence_2(1,i)+error_tmp;
            end
            lms_system(0,'c');
            for i=11:sequence_length
               [~,error_tmp]=lms_system(i);
               error_square_sequence_2(2,i)=error_square_sequence_2(2,i)+error_tmp;
            end            
        end
        error_square_sequence_2=error_square_sequence_2/20;
        figure;
        plot(error_square_sequence_2(1,:),'b');
        hold on;
        plot(error_square_sequence_2(2,:),'r');        
        legend('RLS','LMS')
        title('RLS与LMS算法误差平方收敛对比')
        xlabel('数字时间')
        ylabel('误差平方')
        saveas(gcf,'compare.eps','epsc')        
    end
end