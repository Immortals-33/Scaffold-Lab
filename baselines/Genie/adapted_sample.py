import os
import sys
import torch
import argparse
import numpy as np
import time
import logging
from tqdm import tqdm, trange

from genie.utils.model_io import load_model
from genie.diffusion.genie import Genie
from genie.config import Config

logger = logging.getLogger(__name__)


def main(args):

	# device
	device = 'cuda:{}'.format(args.gpu) if args.gpu is not None else 'cpu'

	# model
	#model = load_model(args.rootdir, args.model_name, args.model_version, args.model_epoch).to(device)
	cfg = Config('weights/scope_l_256/configuration')
	model = Genie.load_from_checkpoint('weights/scope_l_256/epoch=29999.ckpt', config=cfg).to(device)

	# output directory
	#outdir = os.path.join(model.rootdir, model.name, 'version_{}'.format(model.version), 'samples')
	outdir = 'opt'
	#if not os.path.exists(outdir):
	#	os.mkdir(outdir)
	#outdir = os.path.join(outdir, 'epoch_{}'.format(model.epoch))
	#if os.path.exists(outdir):
	#	print('Samples existed!')
	#	sys.exit(0)
	#else:
	#	os.mkdir(outdir)

	# sanity check
	min_length = 50
	max_length = 250
	interval = 50
	max_n_res = model.config.io['max_n_res']
	assert max_length <= max_n_res

	# sample
	for length in trange(min_length, max_length + 1, interval):
		start = time.time()
		for batch_idx in range(args.num_batches):
			mask = torch.cat([
				torch.ones((args.batch_size, length)),
				torch.zeros((args.batch_size, max_n_res - length))
			], dim=1).to(device)
			ts = model.p_sample_loop(mask, args.noise_scale, verbose=False)[-1]
			for batch_sample_idx in range(ts.shape[0]):
				sample_idx = batch_idx * args.batch_size + batch_sample_idx
				coords = ts[batch_sample_idx].trans.detach().cpu().numpy()
				coords = coords[:length]
				np.savetxt(os.path.join(outdir, f'{length}_{sample_idx}.npy'), coords, fmt='%.3f', delimiter=',')
		end = time.time()
		elpased = end - start
		logger.info(f'Time usage for length {length} is {elpased:.3f} and for {elpased / (args.num_batches * args.batch_size):.3f} per backbone.')
		


if __name__ == '__main__':

	# parse arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--gpu', type=str, help='GPU device to use')
	parser.add_argument('-r', '--rootdir', type=str, help='Root directory (default to runs)', default='runs')
	parser.add_argument('-n', '--model_name', type=str, help='Name of Genie model', required=True)
	parser.add_argument('-v', '--model_version', type=int, help='Version of Genie model')
	parser.add_argument('-e', '--model_epoch', type=int, help='Epoch Genie model checkpointed')
	parser.add_argument('--batch_size', type=int, help='Batch size', default=5)
	parser.add_argument('--num_batches', type=int, help='Number of batches', default=2)
	parser.add_argument('--noise_scale', type=float, help='Sampling noise scale', default=0.6)
	args = parser.parse_args()

	# run
	main(args)