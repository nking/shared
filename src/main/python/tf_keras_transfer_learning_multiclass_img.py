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

print(f'tf version={tf.__version__}')

NUM_EPOCHS = 25
BATCH_SZ = 32 # keep between 2 and 32, inclusive
IMG_HEIGHT = 32
IMG_WIDTH = 32
L2_REG = 1E-3
seed = 1234

N_DATALOAD_WORKERS = 0 # defaults to using main thread. uses less memory in total says docs

add_data_augmentation = True

run_small_dataset = True

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

def load_dataset_run_model(plot_title: str ='', useRegularization: bool = False, useDropout: bool = False,
                           add_dense_layers: bool = False, unfreeze_model: bool = False) :

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

    # run the model
    base_model = keras.applications.ResNet50(
            weights='imagenet',  # Load weights pre-trained on ImageNet.
            input_shape=(IMG_WIDTH, IMG_HEIGHT, 3),
            include_top=False, pooling='avg')
    #TODO: consider 'unfreezing' only the last couple of layers
    #TODO: consider ensembling
    if unfreeze_model:
        base_model.trainable = True
    else:
        base_model.trainable = False

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
    x = base_model(inputs, training=False)
    # base_model last layer is 512
    if add_dense_layers:
        # this does not change results much for small subset.  shows should probably unfreeze the pre-trained model
        outputs = layers.Dropout(0.5)(x)
        outputs = layers.Dense(256, activation='elu')(x)
        outputs = layers.Dropout(0.5)(x)
        outputs = layers.Dense(128, activation='elu')(x)
        #outputs = layers.Dropout(0.5)(x)
        #outputs = layers.Dense(64, activation='elu')(x)
    elif useDropout:
        outputs = keras.layers.Dropout(0.5)(x)
    if useRegularization:
        outputs = keras.layers.Dense(num_classes, activation='softmax', kernel_regularizer = keras.regularizers.l2(L2_REG))(x)
    else:
        outputs = keras.layers.Dense(num_classes, activation='softmax')(x)
    model = keras.Model(inputs, outputs)

    model.compile(optimizer='adam',
                  loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=False),
                  metrics=['accuracy'])

    history_train = model.fit(train_ds, validation_data=val_ds,epochs=NUM_EPOCHS)

    plot_util.plot_loss_acc_tf(plot_title+' cifar10', history=history_train)

    loss, acc = model.evaluate(test_ds, verbose=2)
    print(f'{plot_title} cifar10 test: loss={loss:5.2f} accuracy={acc:5.2f}')

if run_small_dataset:
    #load_dataset_run_model("small subset", useRegularization=False)
    #print(f'That was an example of over-fitting')

    #TODO: reduce the over-fitting of the small dataset
    add_data_augmentation = True
    load_dataset_run_model("small subset", useRegularization=True, useDropout=True, add_dense_layers=True,
                           unfreeze_model=True)
    # unfreeze, w/ L2 reg, + dropout is 14 sec/epoch
    #    same + dense layers is

data_dir = data_dir_0
add_data_augmentation = True
load_dataset_run_model("full dataset")

