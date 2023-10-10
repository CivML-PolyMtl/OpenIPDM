% TAGI data processing parameters
split_ratio = 0.85;
shuffle = true;
seed = false;
% misc parameters
reinit_wb = false;
training_plot_save_freq = 1; % used in the 'mod' operator

if ~AnnModel.cont_training
    % TAGI parameters for initializations
    AnnModel.batch_size = 16;
    AnnModel.max_epoches = 30;
    mw_gain = 1;
    Sw_gain = 1;
    mb_gain = 1;
    Sb_gain = 1;
    h_nodes = 128;
    h_layer_act = 4;
    o_layer_act = 0;
    % Initialize the network
    x_dim = size(x_train_val_cont, 2);
    if size(x_train_val_cont, 1) < 200
        AnnModel.batch_size = 2;
    end
    y_dim = 1;
    net = TAGI_init(modelName, AnnModel.max_epoches, AnnModel.batch_size, x_dim, y_dim,...
                    mw_gain, Sw_gain, mb_gain, Sb_gain,...
                    h_nodes, h_layer_act, o_layer_act);
    AnnModel_loaded = [];
else
    AnnModel_loaded = AnnModel;
    net = AnnModel.net;
end