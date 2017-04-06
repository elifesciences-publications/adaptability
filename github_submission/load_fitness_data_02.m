%%%% THIS FILE LOADS DATA RELATED TO THE FITNESS ASSAY PEFORMED BY
%%%% ELIZABETH JERISON IN MAY 2016 (RAW DATA FILES CREATED ON JAN 10, 2017)
%%%% THIS FILE CREATES DATA STRUCTURE KRUG_DATA

krug_data.fnd_hd = {'INIT FIT YPD30','ERR','INIT FIT SC37','ERR'};
krug_data.n_rep = 4;

%% Loading control replicate measurements
filename = '/Users/skryazhi/Dropbox WORK/Dropbox/kruglyak_assays_and_qtl_detection/_SKANALYSIS/data/CSV/control_replicate_measurements.csv';
fid = fopen(filename);
data = textscan(fid, '%s', 'delimiter', '\n');
data = data{1};
fclose(fid);

% tmp = textscan(data{1}, '%s', 'delimiter', ';');
% tmp = tmp{1};

krug_data.n_errest = length(data) - 1;

for ipop = 1:krug_data.n_errest
    tmp = textscan(data{ipop+1}, '%s', 'delimiter', ';');
    tmp = tmp{1}';
    tmp{2} = textscan(tmp{2}, '%s', 'delimiter', ',');
    if ipop == 1
        krug_data.errest = nan( krug_data.n_errest, length(tmp{2}{1}), 2);
    end
    krug_data.errest(ipop,1:length(tmp{2}{1}),1) = str2double(tmp{2}{1})';
    tmp{3} = textscan(tmp{3}, '%s', 'delimiter', ',');
    krug_data.errest(ipop,1:length(tmp{3}{1}),2) = str2double(tmp{3}{1})';
end

clear filename fid data tmp ipop;

%% Loading fitness measurements
% filename = '/Users/skryazhi/Dropbox WORK/Dropbox/kruglyak_assays_and_qtl_detection/_SKANALYSIS/data/CSV/fitness_measurements.csv';
filename = '/Users/skryazhi/Dropbox WORK/Dropbox/kruglyak_adaptability_paper/data 2017-01-10/fitness_measurements_with_population_names_12_29_2016.csv'
fid = fopen(filename);
data = textscan(fid, '%s', 'delimiter', '\n');
data = data{1};
fclose(fid);

tmp = textscan(data{1}, '%s', 'delimiter', ';');
tmp = tmp{1};

krug_data.n_fnd = length(data) - 1;
krug_data.fnd_id = cell( krug_data.n_fnd, 1);
krug_data.fnd_fit = nan( krug_data.n_fnd , 4);
% col 1 and 2 -- mean fit and std in YPD
% col 3 and 4 -- mean fit and std in SC37

krug_data.pop_fit = nan( krug_data.n_fnd , 16);
% col 1 through 8 -- evolved in YPD30; col 9 through 16 -- evolved in SC37;
% col 1:2:7 fit in YPD30; col 2:2:8 fit in SC37
% col 9:2:15 fit in YPD30; col 10:2:16 fit in SC37

for ifnd = 1:krug_data.n_fnd
    tmp = textscan(data{ifnd+1}, '%s', 'delimiter', ';');
    tmp = tmp{1};
    krug_data.fnd_id{ifnd} = tmp{1};
    
    krug_data.fnd_fit(ifnd,1:4) = str2double( tmp(2:5) );
    
    tmp{6} = textscan( tmp{6}, '%s', 'delimiter', ',');    
    for ipop = 1:length(tmp{6}{1})
        tmp1 = textscan( tmp{6}{1}{ipop}, '%s', 'delimiter', ' ');
        tmp1 = tmp1{1};
        tmp2 = regexp(tmp1{1}, 'clone-(\d)', 'tokens');
        pop_id = str2double(tmp2{1});
        k = 2 * (pop_id-1) + 1;
        krug_data.pop_fit(ifnd, k:k+1) = str2double( tmp1(2:3) );
    end
    
    tmp{7} = textscan( tmp{7}, '%s', 'delimiter', ',');    
    for ipop = 1:length(tmp{7}{1})
        tmp1 = textscan( tmp{7}{1}{ipop}, '%s', 'delimiter', ' ');
        tmp1 = tmp1{1};
        tmp2 = regexp(tmp1{1}, 'clone-(\d)', 'tokens');
        pop_id = str2double(tmp2{1});
        k = 2 * (pop_id-1) + 9;
        krug_data.pop_fit(ifnd, k:k+1) = str2double( tmp1(2:3) );
    end
    

end

clear filename fid data tmp ipop ifnd tmp1 k tmp2 pop_id;



%% Loading founder genotypes
% filename = '/Users/skryazhi/Dropbox WORK/Dropbox/kruglyak_assays_and_qtl_detection/_SKANALYSIS/data/CSV/segregant_genotypes.csv';
filename = '/Users/skryazhi/Dropbox WORK/Dropbox/kruglyak_adaptability_paper/data 2017-01-10/segregant_genotypes_with_header.csv';
fid = fopen(filename);
data = textscan(fid, '%s', 'delimiter', '\n');
data = data{1};
fclose(fid);

tmp = textscan(data{1}, '%s', 'delimiter', ',');
tmp = tmp{1};
krug_data.snp_id = tmp(2:end);
krug_data.n_snp = length( krug_data.snp_id );



krug_data.fnd_gt =  false( krug_data.n_fnd, krug_data.n_snp);

for ifnd = 2:krug_data.n_fnd
    tmp = textscan(data{ifnd}, '%s', 'delimiter', ';');
    tmp = tmp{1};
    
    ix = find( strcmp( tmp{1}, krug_data.fnd_id ) );
        
    tmp{2} = textscan( tmp{2}, '%s', 'delimiter', ',');
    tmp1 = str2double( tmp{2}{1} );
    krug_data.fnd_gt(ix, :) = logical( tmp1' );
end

clear filename fid data tmp ipop ifnd tmp1 ix;

% filename = '/Users/skryazhi/Dropbox WORK/Dropbox/kruglyak_assays_and_qtl_detection/_SKANALYSIS/data/krug_data_01.mat';
filename = '/Users/skryazhi/Dropbox WORK/Dropbox/kruglyak_adaptability_paper/data 2017-01-10/krug_data_02.mat';
save(filename, 'krug_data');