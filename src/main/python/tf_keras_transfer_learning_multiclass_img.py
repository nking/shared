# multiclass classification through transfer learning as feature extractor
# script uses tenserflow and keras
# dataset=cifar10
# transfer model=ResNet50 pretrained on imagenet

# adapted from tensorflow.org  and keras tutorials

# Author: Nichole L. King

import tempfile
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from PIL import Image
from tempfile import TemporaryDirectory
import tensorflow as tf
from tensorflow import keras
import plot_util
import tf_dataset_util
from keras import layers
from tensorflow.math import confusion_matrix

print(f'tf version={tf.__version__}')

NUM_EPOCHS = 25
BATCH_SZ = 32 # keep between 2 and 32, inclusive
IMG_HEIGHT = 32
IMG_WIDTH = 32
L2_REG = 1E-3 # 1E-3 to 0.1
DO_PROB = 0.5 # 0 to 1.  drop out probability.
seed = 1234

#TODO: consider ensembling when accuracy is high enough that a few percent increase matters


N_DATALOAD_WORKERS = 0 # defaults to using main thread. uses less memory in total says docs

add_data_augmentation = True

run_small_dataset = False

##TODO: replace / with os.path.sep

if True:
    # desktop environment:
    data_dir = os.environ["ML_DATASETS_HOME"] + '/cifar10png'
else:
    # colab environment
    # other args: url, folder_in_archive
    data_dir = tempfile.mkdtemp()
    print(f'temp_dataset_dir={data_dir}')
    #!pip install cifar2png
    #%cifar2png cifar10 data_dir

if run_small_dataset:
    # make a small subset of the data.  copies n_test = n_train by default

    add_data_augmentation = False

    # Error with this one.  TODO: fix shape
    # ValueError: Input 0 of layer "model" is incompatible with the layer: expected shape=(None, 32, 32, 3),
    # found shape=(32, 32, 3)
    #[train_ds, val_ds, class_names] = tf_dataset_util.get_subset_train_val_datasets_0(image_size=(IMG_WIDTH, IMG_HEIGHT),
    #                               num_channels=3, crop_to_aspect_ratio=False, batch_size=BATCH_SZ,
    #                              data_dir=data_dir, n_train=120, fraction_val=0.8, seed=seed, verbose=1)

    data_dir_0 = data_dir
    data_dir_small = tf_dataset_util.copy_subset_to_TMP(data_dir=data_dir, n_train=120, n_test=120, fraction_val=0.8, seed=seed, verbose=1)
    data_dir = data_dir_small

def load_dataset_run_model(plot_title: str ='', use_regularization: bool = False, use_dropout: bool = False,
                           n_unfreeze:int=0, unfreeze_all:bool=False, add_dense_layers:bool=False) :

    # type tf.data.Dataset
    train_ds, val_ds = tf.keras.utils.image_dataset_from_directory(
        data_dir + '/train', labels="inferred",
        label_mode="int", class_names=None,
        color_mode="rgb", batch_size=BATCH_SZ,image_size=(IMG_HEIGHT, IMG_WIDTH),
        shuffle=True, seed=seed,
        validation_split=0.2,
        subset="both",
        interpolation="bilinear", follow_links=False, crop_to_aspect_ratio=False,
    )
    test_ds = tf.keras.utils.image_dataset_from_directory(
        data_dir + '/test', labels="inferred",
        label_mode="int", class_names=None,
        color_mode="rgb", batch_size=BATCH_SZ,image_size=(IMG_HEIGHT, IMG_WIDTH),
        shuffle=True, seed=seed,
        validation_split=None,
        interpolation="bilinear", follow_links=False, crop_to_aspect_ratio=False,
    )
    train_path = data_dir + '/train'
    #num_classes = sum([1 for entry in os.scandir(train_path) if not entry.name.startswith('.') and entry.is_dir()])
    num_classes = len(set(train_ds.class_names))
    print(f'num_classes={num_classes}')

    #plot_util.plot_first_tf(train_ds, train_ds.class_names)
    plot_util.print_ds_shape_tf(train_ds)

    # build, compile, run the model

    # include_Top = False does not add the last layer to the model
    #     last layer would have been
    #     x = layers.Dense(classes, activation=classifier_activation, name="predictions")(x)
    # we supply that below, tailored to train_ds
    base_model = keras.applications.ResNet50(
            weights='imagenet',  # Load weights pre-trained on ImageNet.
            input_shape=(IMG_WIDTH, IMG_HEIGHT, 3),
            include_top=False, pooling='avg')
            # if include_top=True, add , classifier_activation='softmax'

    #TODO: consider ensembling

    base_model.trainable = False
    if unfreeze_all:
        base_model.trainable = True
    else:
        print(f'base model {base_model.name} has {len(base_model.layers)} layers')
        # base model has 176 layers
        if n_unfreeze > len(base_model.layers):
            raise ValueError(f"model has {len(base_model.layers)}.  too large n_unfreeze={n_unfreeze}")
        for i in range(n_unfreeze):
            base_model.layers[-(i+1)].trainable = True

    if add_data_augmentation:
        data_augmentation = keras.Sequential(
            [layers.RandomFlip("horizontal", input_shape=(IMG_WIDTH, IMG_HEIGHT, 3)),
              layers.RandomRotation(0.1),
              layers.RandomZoom(0.1),]
        )

    inputs = keras.Input(shape=(IMG_WIDTH, IMG_HEIGHT, 3))
    inputs = keras.applications.resnet50.preprocess_input(inputs)
    if add_data_augmentation:
        inputs = data_augmentation(inputs)
    x = base_model(inputs)
    # base_model last layer is 512

    if add_dense_layers:
        # this does not change results much for small subset.  shows should probably unfreeze the pre-trained model
        outputs = layers.Dropout(DO_PROB)(x)
        outputs = layers.Dense(256, activation='elu')(x)
        outputs = layers.Dropout(DO_PROB)(x)
        outputs = layers.Dense(128, activation='elu')(x)
        outputs = layers.Dropout(DO_PROB)(x)
        outputs = layers.Dense(64, activation='elu')(x)
        outputs = layers.Dropout(DO_PROB)(x)
        outputs = layers.Dense(DO_PROB, activation='elu')(x)
    elif use_dropout:
        outputs = keras.layers.Dropout(DO_PROB)(x)

    if use_regularization:
        outputs = keras.layers.Dense(num_classes, activation='softmax',
                                     kernel_regularizer = keras.regularizers.l2(L2_REG), name='predictions')(x)
    else:
        outputs = keras.layers.Dense(num_classes, activation='softmax', name='predictions')(x)

    model = keras.Model(inputs, outputs)

    model.compile(optimizer='adam',
                  loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=False),
                  metrics=['accuracy'])

    history_train = model.fit(train_ds, validation_data=val_ds,epochs=NUM_EPOCHS)

    loss, acc = model.evaluate(test_ds, verbose=0)
    print(f'{plot_title} cifar10 test: loss={loss:5.2f} accuracy={acc:5.2f}')

    plot_util.plot_loss_acc_tf(plot_title+' cifar10', history=history_train, test_acc=acc, test_loss=loss)

    # ===== calc precision, recall, and F1  ======
    predictions = model.predict(test_ds)
    #print(f'predictions.shape={predictions.shape}')
    # shape is nData x nClasses
    n_test_data = predictions.shape[0]
    l0 = np.zeros((n_test_data), dtype="int32")
    p0 = np.zeros((n_test_data), dtype="int32")
    j = 0
    for image_batch, labels_batch in test_ds:
        for i in range(image_batch.shape[0]):
            l0[j] = labels_batch[i]
            p0[j] = np.argmax(predictions[j])
            j += 1
    # The confusion matrix columns represent the prediction labels and the rows represent the
    #    real labels.
    cm = confusion_matrix(l0, p0)

    # https://dev.to/overrideveloper/understanding-the-confusion-matrix-264i
    rowsums = np.sum(cm, axis=0)
    colsums = np.sum(cm, axis=1)
    cmsum = sum(rowsums)
    diagsum = np.trace(cm) # scalar
    # tp[i] = cm[i][i]
    # tn[i] = cmsum - rowsums[i] - solsums[i]
    # fp[i] = colsums[i] - cm[i][i]
    # fn[i] = rowsums[i] - cm[i][i]

    # accuracy = diagsum / cmsum
    # precision[i] = tp[i]/(tp[i] + fp[i]) = np.trace(cm) / colsums
    # recall[i] = tp[i]/(tp[i] + fn[i])    = np.trace(cm) / rowsums
    # f1[i] = 2./( (1./precision[i]) + (1./recall[i])) = ...
    # macro_precision = sum(precision) / (len(precision))
    # macro_recall = sum(recall) / (len(recall))

    micro_precision = diagsum / sum(colsums)
    micro_recall = diagsum / sum(rowsums)
    micro_f1 = 2./((1./micro_precision) + (1./micro_recall))
    accuracy2 = diagsum/cmsum
    # ideally, prec and rec > 0.8 and F1 close to 1
    print(f'accuracy2={accuracy2:.3f}, precision={micro_precision:.3f}, recall={micro_recall:.3f}, f1={micro_f1:.3f}')

if False and run_small_dataset:
    load_dataset_run_model("small subset (over-fitting)", use_regularization=False)
    print(f'That was an example of over-fitting')

    add_data_augmentation = True
    # retrain all layers runtime is 20-40 sec/epoch
    load_dataset_run_model(f'small (n_unfreeze=0, DropOutP={DO_PROB}, DataAug={add_data_augmentation})', n_unfreeze=0, use_dropout=True)
    load_dataset_run_model(f'small (n_unfreeze=50, DropOutP={DO_PROB}, DataAug={add_data_augmentation})', n_unfreeze=50, use_dropout=True)
    load_dataset_run_model(f'small (unfreeze all, DropOutP={DO_PROB}, DataAug={add_data_augmentation})', unfreeze_all=True, use_dropout=True)

if not run_small_dataset:
    data_dir = data_dir_0
    add_data_augmentation = True
    dlabel = "full"
    load_dataset_run_model(f'{dlabel} (n_unfreeze=0, DropOutP={DO_PROB}, DataAug={add_data_augmentation})', use_dropout=True, n_unfreeze=0)
    print(f'a model that is not over-fitting was just plotted')
else:
    dlabel = "small"
    NUM_EPOCHS = 5
    L2_REG = 1E-2
    load_dataset_run_model(f'{dlabel} (n_unfreeze=0, DropOutP={DO_PROB}, L2_REG={L2_REG:.1e}, DataAug={add_data_augmentation})',
                           use_dropout=True, n_unfreeze=0, use_regularization=True)
    L2_REG = 1E-1
    load_dataset_run_model(f'{dlabel} (n_unfreeze=0, DropOutP={DO_PROB}, L2_REG={L2_REG:.1e}, DataAug={add_data_augmentation})',
                           use_dropout=True, n_unfreeze=0, use_regularization=True)

