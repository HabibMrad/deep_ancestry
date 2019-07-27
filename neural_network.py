import tensorflow as tf
tf.logging.set_verbosity(tf.logging.ERROR)

import numpy as np
import gzip
from argparse import ArgumentParser


# DEFAULT: ReLU(inner), Softmax(out), cross-entropy, Adam
# parse arguments
parser = ArgumentParser(description='Mandatory commands: file path and at least a population name')
parser.add_argument('-file', type=str, help='Path and name prefix of the infiles, whose name format is assumed to be [path + name_prefix](_pops).(train/test).gz')
parser.add_argument('-epoch', type=int, default=1, help='Number of epochs to run on the same training set')
parser.add_argument('-inner', type=int, default=0, help='Number of inner layers')
parser.add_argument('-neurons', type=int, nargs='*', default=[1000], help='List of the number of neurons per inner layer, ordered')
parser.add_argument('-drop', type=float, default=0., help='Proportion of neurons to drop out in each training run')
parser.add_argument('-trbatch', type=int, default=50, help='Number of simulated samples per training batch')
parser.add_argument('-tsbatch', type=int, default=10000, help='Number of simulated samples per testing batch')
parser.add_argument('-optim', type=str, default='Adam', help='Optimizer algorithm, Adam by default, GradientDescent optional')
parser.add_argument('-rate', type=float, default=0.001, help='Learning rate')
parser.add_argument('-cost', type=float, default=0., help='Cost that, if reached, causes the training to finish')
parser.add_argument('-cnst', type=float, default=1e-10, help='Tiny constant added to limit low extremes in soft label of admixed individuals (default 1e-10) so there are no log(0) == Nan. Should be the same in offspring_simulator and in neural_network')
parser.add_argument('-display', type=int, default=10, help='Batch interval to show')
parser.add_argument('-log', type=str, default=False, help='Path for log files')
parser.add_argument('-seed', type=int, default=False, help='Set random seed')
args = parser.parse_args()

file = args.file

if args.seed:
    tf.set_random_seed(args.seed)

with gzip.open(f"{file}.train.gz", "rt") as training_genotypes, \
    gzip.open(f"{file}_pops.train.gz", "rt") as training_pops:

    # Set parameters
    epochs = args.epoch

    neurons_input = int(next(training_genotypes))
    num_inner = args.inner
    neurons_inner = args.neurons
    neurons_output = int(next(training_genotypes))

    dropout = args.drop

    num_simulations = int(next(training_genotypes))
    training_batch_size = args.trbatch
    num_training_batches = int(num_simulations / training_batch_size)
    optim = args.optim
    learning_rate = args.rate
    display_batch = args.display
    testing_batch_size = args.tsbatch
    final_cost = args.cost

    # TF graph input
    x = tf.placeholder("float", [None, neurons_input])    # number input neurons (= num SNPs with genotype info provided as input)
    y = tf.placeholder("float", [None, neurons_output])   # number output neurons (= num pops)
    z_dropout = tf.placeholder(tf.float32)

    ## Create a model
    # Set model weights and bias neuron for output layer
    if num_inner == 0:
        # there is only an out layer
        W_out = tf.get_variable("weights_out", shape=[neurons_input, neurons_output], initializer=tf.glorot_uniform_initializer())
        b_out = tf.get_variable("bias_out", shape=[neurons_output], initializer=tf.constant_initializer(0.1))

        # Construct a linear model for output layer
        with tf.name_scope("W_out_x_b") as scope:
            if dropout == 0.:
                model_out = tf.nn.softmax(tf.matmul(x, W_out) + b_out)  # Softmax
            else:
                # Dropout
                drop = tf.nn.dropout(x, rate=z_dropout)
                model_out = tf.nn.softmax(tf.matmul(drop, W_out) + b_out) # Softmax
    else:
        # there are inner layers
        layers = dict()
        for layer in range(1, num_inner + 1):
            ## inner layers biases, weights and models
            # biases
            layers[f"b{layer}"] = tf.get_variable(f"bias{layer}", shape=[neurons_inner[layer-1]], initializer=tf.constant_initializer(0.1))
            # weigths and models
            if layer == 1:
                # 1st dimension of W is the number of input neurons
                layers[f"W{layer}"] = tf.get_variable(f"weights{layer}", shape=[neurons_input, neurons_inner[layer - 1]], initializer=tf.glorot_uniform_initializer())
                # Construct a logistic model for inner layers
                with tf.name_scope(f"W{layer}_x_b") as scope:
                    if dropout == 0.:
                        layers[f"model{layer}"] =  tf.nn.leaky_relu(tf.matmul(x, layers[f"W{layer}"]) + layers[f"b{layer}"])  # leaky ReLU
                    else:
                        # Dropout
                        drop = tf.nn.dropout(x, rate=z_dropout)
                        layers[f"model{layer}"] =  tf.nn.leaky_relu(tf.matmul(drop, layers[f"W{layer}"]) + layers[f"b{layer}"])  # leaky ReLU
            else:
                # 2nd to last inner layers
                # 1st dimension of W is the number of neurons of the previous layer
                layers[f"W{layer}"] = tf.get_variable(f"weights{layer}", shape=[neurons_inner[layer - 2], neurons_inner[layer - 1]], initializer=tf.glorot_uniform_initializer())
                # Construct a logistic model for inner layers
                with tf.name_scope(f"W{layer}_x_b") as scope:
                    if dropout == 0.:
                        layers[f"model{layer}"] =  tf.nn.leaky_relu(tf.matmul(layers[f"model{layer - 1}"], layers[f"W{layer}"]) + layers[f"b{layer}"])  # leaky ReLU
                    else:
                        # Dropout
                        drop = tf.nn.dropout(layers[f"model{layer - 1}"], rate=z_dropout)
                        layers[f"model{layer}"] =  tf.nn.leaky_relu(tf.matmul(drop, layers[f"W{layer}"]) + layers[f"b{layer}"])  # leaky ReLU
            if args.log:
                # Add summary ops to collect data
                layers[f"w{layer}_h"] = tf.summary.histogram(f"weights{layer}", layers[f"W{layer}"])
                layers[f"b{layer}_h"] = tf.summary.histogram(f"biases{layer}", layers[f"b{layer}"])

        # rename the last model for calling it unambiguously in the output layer
        model_last_inner = layers[f"model{layer}"]
        del layers[f"model{layer}"]

        # extract all weights, biases and models as independent variables
        locals().update(layers)

        # output layer
        W_out = tf.get_variable("weights_out", shape=[neurons_inner[-1], neurons_output], initializer=tf.glorot_uniform_initializer())
        b_out = tf.get_variable("bias_out", shape=[neurons_output], initializer=tf.constant_initializer(0.1))
        # Construct a linear model for output layer
        with tf.name_scope("W_out_x_b") as scope:
            if dropout == 0.:
                model_out = tf.nn.softmax(tf.matmul(model_last_inner, W_out) + b_out) # Softmax
            else:
                # Dropout
                drop = tf.nn.dropout(model_last_inner, rate=z_dropout)
                model_out = tf.nn.softmax(tf.matmul(drop, W_out) + b_out) # Softmax

    if args.log:
        # Add summary ops to collect data from output layer
        w_out_h = tf.summary.histogram("weights_out", W_out)
        b_out_h = tf.summary.histogram("biases_out", b_out)

    # calculate loss (cost function)
    with tf.name_scope("cost_function") as scope:
        # categorical cross-entropy for one-hot labels, i.e. non-admixed individuals
        # note: it works also in case true labels are soft probabilities (e.g. [0.8, 0.2, 0] = 20% CEU 80% YRI admixed ind), although the minimum loss does not reach 0,
        #       and this is higher the more spread are the values e.g. loss(y=ŷ=[0.5, 0.5, 0.5]) > loss(y=ŷ=[1, 0, 0]) = 0,
        #       but with batches large enough this should be normalized
        # note: this outputs the AVERAGE loss across all simulations within a batch, so the batch size does not change the scale of the loss
        #       (which could slow things a bit)
        # note: tiny constant added to limit low extremes so they do not reach 0, since log(0) == Nan
        cnst = args.cnst
        cost_function = tf.reduce_mean(abs(tf.reduce_sum(y * tf.log(model_out + cnst), axis=1)))
        # note: the thing is to train using pure inds and then test with admixed,
        # using the same loss function. If this does not work, maybe train with admix indivs?

        if args.log:
            # Create a summary to monitor the cost function
            tf.summary.scalar("cost_function", cost_function)

    # Minimize loss using optimizer (performs gradient descent + back-propagation)
    with tf.name_scope("train") as scope:
        # Adam optimizer
        if optim == 'Adam':
            optimizer = tf.train.AdamOptimizer(learning_rate).minimize(cost_function)
        elif optim == 'GradientDescent':
            # regular GradientDescent
            optimizer = tf.train.GradientDescentOptimizer(learning_rate).minimize(cost_function)

    # Initializing the variables
    init = tf.global_variables_initializer()

    if args.log:
        # Merge all summaries into a single operator
        merged_summary_op = tf.summary.merge_all()


    ## Launch the graph
    with tf.Session() as sess:
        sess.run(init)

        if args.log:
            # save logs
            summary_writer = tf.summary.FileWriter(args.log, graph=sess.graph)

        # print log message
        print(f"\n\nRunning {epochs} epochs using {neurons_input} SNPs, {num_inner} inner layers of {neurons_inner} neurons each, dropping out {dropout * 100}% neurons, guessing {neurons_output} populations, applying a learning rate of {learning_rate} on a batch size of {training_batch_size} over {num_simulations} training simulations until reaching a cost = {final_cost}, using the {optim}Optimizer algorithm. Then, testing using a batch size of {testing_batch_size} simulations...\n")

        # Training cycle (nº trained ANN)
        for iteration in range(epochs):
            if iteration > 0:
                training_genotypes.seek(0)
                [next(training_genotypes) for i in range(3)]
                training_pops.seek(0)
            epoch_cost = 0.
            # Loop over all batches
            for batch in range(1,num_training_batches + 1):

                # store batch size simulations within a batch
                batch_genos = np.genfromtxt(training_genotypes, max_rows=training_batch_size, dtype=float)
                batch_pops = np.genfromtxt(training_pops, max_rows=training_batch_size, dtype=str)

                # Fit training using batch data
                try:
                    _, cost = sess.run([optimizer, cost_function], feed_dict={x: batch_genos, y: batch_pops, z_dropout: dropout})
                except ValueError:
                    # batches consist of a single simulation
                    batch_genos = [batch_genos]
                    batch_pops = [batch_pops]
                    _, cost = sess.run([optimizer, cost_function], feed_dict={x: batch_genos, y: batch_pops, z_dropout: dropout})
                # Compute the average loss per simulation within the batch
                epoch_cost += cost

                if args.log:
                    # Write logs for each iteration
                    summary_str = sess.run(merged_summary_op, feed_dict={x: batch_genos, y: batch_pops, z_dropout: dropout})
                    summary_writer.add_summary(summary_str, iteration * num_training_batches + batch)

                # Display logs per batch
                if batch % display_batch == 0 or batch == 1:
                    print(f"Batch nº{batch} avg.cost = {'%.9f' % cost}; epoch current cost = {epoch_cost / batch}\n")

                # finish training if reached error
                if cost <= final_cost:
                    break
            # Display logs per ANN
            print(f"Finished epoch nº{iteration + 1}; final batch avg.cost = {'%.9f' % cost}; epoch cost = {epoch_cost / batch}\n")

        ## Test the ANN and calculate accuracy
        print(f"Testing accuracy...\n")

        predictions = tf.equal(tf.argmax(model_out, 1), tf.argmax(y, 1))
        accuracy = tf.reduce_mean(tf.cast(predictions, "float"))

        with gzip.open(f"{file}.test.gz", "rt") as testing_genotypes, \
            gzip.open(f"{file}_pops.test.gz", "rt") as testing_pops:

            batch_genos = np.genfromtxt(testing_genotypes, max_rows=testing_batch_size, dtype=float)
            batch_pops = np.genfromtxt(testing_pops, max_rows=testing_batch_size, dtype=str)

            print("Accuracy:", accuracy.eval({x: batch_genos, y: batch_pops, z_dropout: 0.}))
