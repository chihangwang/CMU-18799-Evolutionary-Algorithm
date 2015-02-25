"""
	Student Name: Chih-Ang Wang
	AndrewID: chihangw
	README:  This python performs evolutionary strategy on Ackley function.
			 The best (minimal) value found with seed=1234 is 5.1639e-08, which
			 can be verified by running this program (python es.py) in your
			 terminal. The output string is in the form of:
			 [best_ackley_val, [x1, x2,..., x30, sigma1, sigma2,... sigma30]]
"""
import random
import sys
import math

class EvoStrategy:

	def __init__(self, seed, n):
		random.seed(seed)
		self.n = n
		self.epsilon = 0.01
		self.tau = (1/math.sqrt(2*n))*2
		self.pools = [[]]*30
		self.offsprings = [[]]*200
		self.rank_table = []

	def init(self):
		for i in range(30):
			# 30 object variables, 30 strategy parameters
			# chromosome = [0] * 60
			chromosome = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,]

			for j in range(30):
				chromosome[j] = round(random.uniform(-30, 30), 8)
			for j in range(30, 60):
				chromosome[j] = 1

			self.pools[i] = chromosome

	def recombine(self):
		# object variables(30) ----> local discrete recombination
		# strategy parameters(30) -> global intermediary recombination
		for cnt in range(200):
			index1 = int(math.floor(random.random()*30))
			index2 = int(math.floor(random.random()*30))
			parent1 = self.pools[index1]
			parent2 = self.pools[index2]

			child = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,]

			# local discrete recombination
			for i in range(30):
				dice = random.getrandbits(1)
				if dice == 0:
					child[i] = round(parent1[i], 8)
				else:
					child[i] = round(parent2[i], 8)

			# global intermediary recombination
			for i in range(30, 60):
				index1 = int(math.floor(random.random()*30))
				index2 = int(math.floor(random.random()*30))
				parent1 = self.pools[index1]
				parent2 = self.pools[index2]
				child[i] = round(0.5*(parent1[i] + parent2[i]), 8)

			self.offsprings[cnt] = child

	def mutate(self):
		# mutate the strategy parameters (sigma) first based on:
		# sigma'(i) = sigma(i) * exp(tau * Ni(0,1))
		for chromosome in self.offsprings:
			N = random.gauss(0, 1)
			for i in range(30, 60):
				Ni = random.gauss(0, 1)
				sigma = chromosome[i] * math.exp(self.tau*N)
				# if sigma < self.epsilon:
				# 	sigma = self.epsilon
				chromosome[i] = round(sigma, 8)

				# then mutate the object variables based on:
				# x'(i) = x(i) + sigma(i) * Ni(0, 1)
				x = chromosome[i-30] + chromosome[i]*Ni
				chromosome[i-30] = round(x, 8)

	def eval_ackley(self):
		# for child in offsprings, evaluate its ackley val and insert back into
		# the offsprings list with [ackley_val, [x1, x2,..., x30, sigma1, ...]]
		for index in range(200):
			child = self.offsprings[index]
			sigma_xi_square = 0
			sigma_cos_2_pi_xi = 0
			for i in range(30):
				xi = child[i]
				sigma_xi_square += math.pow(xi, 2)
				sigma_cos_2_pi_xi += math.cos(2*math.pi*xi)

			val = -20 * math.exp(-0.2*math.sqrt(sigma_xi_square/self.n)) - math.exp(sigma_cos_2_pi_xi/self.n) + 20 + math.e

			pair = [val, child]
			self.offsprings[index] = pair

	def select_survivor(self):
		# rank the offspring, and select the 30 childs with lowest ackley_val
		# as survivors.
		def predicate(child):
			return child[0]

		# sort the offspring list based on ackley_val
		rank_list = sorted(self.offsprings, key=predicate)

		# selection based on (30, 200)
		for i in range(30):
			self.pools[i] = rank_list[i][1]

		# store the best (minimal) value
		self.rank_table.append(rank_list[0])

	def print_best(self):
		print "BEST ===============================================\n"
		def predicate(child):
			return child[0]
		rank_list = sorted(self.rank_table, key=predicate)
		print rank_list[0]


	def show(self):
		print "POOLS =====================================\n"
		for i in range(30):
			print self.pools[i]
		print "\nCHILD =====================================\n"
		for i in range(200):
			print self.offsprings[i]
		print "\nEND =======================================\n\n\n"

def main():
	seed = 1234
	n = 30

	evo = EvoStrategy(seed, n)
	evo.init()

	for i in range(1000):
		evo.recombine()
		evo.mutate()
		evo.eval_ackley()
		evo.select_survivor()

	evo.print_best()
	# evo.show()

if __name__ == '__main__':
	main()