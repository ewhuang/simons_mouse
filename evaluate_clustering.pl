## This script takes clustering output and the network file to compute the in or out density as well as other evaluation measurement for the output.

## It is been tested that it produced the same value as evaluate_clustering-v4.pl in hal

use strict;

die "Usage: perl $0 ClusteringOutput RealNetwork\n" unless @ARGV==2;

my $output = shift (@ARGV); ## one file with given clustering (output from program)
my $realnw = shift(@ARGV); ## one file with real network

my %outclus;
open (OUT, "<$output");
my $hdr = <OUT>;
while (<OUT>) {
    my ($d1, $species, $d2, $gene, $d3, $cluster) = split(/\s+/);
    # warn "$species\t$gene\t$cluster\n";
    $outclus{$species}{$gene} = $cluster; ## outclus{0}{200} = 29
}
close(OUT);

my %network_total_edge_2x;
my %realnetwork;
my %nodeDegree; ## degree of node in a species
my %total_spc_node_degree;
open (REAL, "<$realnw");
$hdr = <REAL>;
while (<REAL>) {
    my ($species, $g1, $g2, $wt) = split(/\s+/);
    $realnetwork{$species}{$g1}{$g2} = $wt; ## realnetwork{0}{559}{365} = 1, this network only counts within species edges, NOT orthologous edges
    $nodeDegree{$species}{$g1} ++;
    $network_total_edge_2x{$species} ++;
}
close(REAL);

for my $spc (keys %nodeDegree){
	for my $g (keys %{$nodeDegree{$spc}}){
		$total_spc_node_degree{$spc} += $nodeDegree{$spc}{$g};
	}
}

print "Descriptive statistics of predicted clusters\n";
DescriptiveStats(\%outclus);

sub DescriptiveStats {
    my ($clus1ptr) = @_;
    my %cs;
    my $clusterNum = 0;
    foreach my $spc (keys %{$clus1ptr}) {
		my @clusters = values %{$$clus1ptr{$spc}};
		foreach my $c (@clusters) { 
			$cs{$c} = 1; ## unique cluster c = 1,2,3 ...
		}
		$clusterNum = keys %cs; ## cluster number
    }
    
    my $sumindensity = 0;
    my $countindensity = 0;
    my $sumoutdensity = 0;
    my $countoutdensity = 0;
    my $suminoutratio = 0;
    my $sum_ave_modularity = 0;
    my %total_cluster_edge = ();
	my %total_cluster_modularity = ();
	my %total_cluster_modularity_weighted = ();
	my %total_cluster_modularity_bySize = ();
	my %total_cluster_degree = ();
	my %cluster_size = (); ## cluster size
	my $sum_modularity_cluster_spc = 0;
	my $sum_modularity_cluster_spc_weighted = 0;
	my $sum_modularity_cluster_spc_bySize = 0;

	my %indensity = ();
	my %outdensity = ();
	my %outdegree = ();
	my %incount = ();
	my %outcount = ();
	my %weighted_in_out_ratio = ();
	my $sum_weighted_in_out_ratio = 0;

	# my $sum_product_verDen_verIntro = 0;

	## for each Cluster
    foreach my $c (keys %cs) {
		## for each Species
		foreach my $spc (keys %{$clus1ptr}) {

		    my @genes = keys %{$$clus1ptr{$spc}};
		    $cluster_size{$spc}{$c} = 0;

		    $indensity{$spc}{$c} = 0;
		    $outdensity{$spc}{$c} = 0;
		    $incount{$spc}{$c} = 0;
		    $outcount{$spc}{$c} = 0;
		   
		    # my $verIntrovert = 0;
		    # my $product_verDen_verIntro = 0;
		    # my $ave_modularity = 0;

		    foreach my $g1 (@genes) {
		    	## g1 in c
				if ($$clus1ptr{$spc}{$g1} eq $c) { 
					## sum genes in a cluster
					$cluster_size{$spc}{$c}++;
					## sum total degree in cluster
			    	$total_cluster_degree{$spc}{$c} += $nodeDegree{$spc}{$g1};

			    	my $inDegree = 0;
			    	foreach my $g2 (@genes) {
			    		## i /= j
						if ($g1 eq $g2) { next; } 
						
						## g2 in c
						if ($$clus1ptr{$spc}{$g2} eq $c) { 

				    		## if there is an edge i,j
				    		if (defined($realnetwork{$spc}{$g1}{$g2})) {
				    			## increase total edges in cluster  
				    			$indensity{$spc}{$c}++;
				    			$inDegree ++;
				    			## sum total edges in cluster
				    			$total_cluster_edge{$spc}{$c} ++;
				    			## sum total modularity in cluster
				  				$total_cluster_modularity{$spc}{$c} += (1 / $network_total_edge_2x{$spc}) - (($nodeDegree{$spc}{$g1} * $nodeDegree{$spc}{$g2}) / ($network_total_edge_2x{$spc} ** 2));
				    		}
				    		else{
				    			## add penalty to nodes without edges
				    			$total_cluster_modularity{$spc}{$c} += - (($nodeDegree{$spc}{$g1} * $nodeDegree{$spc}{$g2}) / ($network_total_edge_2x{$spc} ** 2));
				    		}
				    		## sum total potential edges in-cluster
				    		$incount{$spc}{$c}++; 
						}
						else { ## g2 not in c
				    		if (defined($realnetwork{$spc}{$g1}{$g2})) {
				    			## sum total out-edges
				    			$outdensity{$spc}{$c}++;
				    			$outdegree{$spc}{$c} ++; 
				    		}
				    		## sum total potential edges out-cluster
				    		$outcount{$spc}{$c}++; 
						}
			    	}
	
				}
		    }
		    ## weight modularity by "#total edges in cluster" / "#total degree in cluster"
		    if ($total_cluster_degree{$spc}{$c} >= 1) {
		    	$total_cluster_modularity_bySize{$spc}{$c} = $total_cluster_modularity{$spc}{$c} / $cluster_size{$spc}{$c};
		    	$total_cluster_modularity_weighted{$spc}{$c} = $total_cluster_modularity{$spc}{$c} * ($total_cluster_edge{$spc}{$c} / $total_cluster_degree{$spc}{$c});
		    	
		    }
		    else {
		    	$total_cluster_modularity_bySize{$spc}{$c} = $total_cluster_modularity{$spc}{$c} / 1;
		    	$total_cluster_modularity_weighted{$spc}{$c} = $total_cluster_modularity{$spc}{$c} * ($total_cluster_edge{$spc}{$c} / 1);
		    }
		    ## sum modularity over clusters and spc
		    $sum_modularity_cluster_spc += $total_cluster_modularity{$spc}{$c};
		    $sum_modularity_cluster_spc_bySize += $total_cluster_modularity_bySize{$spc}{$c};
		    $sum_modularity_cluster_spc_weighted += $total_cluster_modularity_weighted{$spc}{$c};

		    ## total nodes in a cluster >= 2
		    if ($incount{$spc}{$c} > 0) {
		    	$indensity{$spc}{$c} = $indensity{$spc}{$c} / $incount{$spc}{$c};
		    	$outdegree{$spc}{$c} = $outdegree{$spc}{$c} / $incount{$spc}{$c};
		    } 

		    if ($outcount{$spc}{$c} > 0) { 
		    	$outdensity{$spc}{$c} = $outdensity{$spc}{$c} / $outcount{$spc}{$c}; 
		    } 

		    if ($incount{$spc}{$c} > 0){
		    	$weighted_in_out_ratio{$spc}{$c} = ($indensity{$spc}{$c} / ($indensity{$spc}{$c} + $outdensity{$spc}{$c})) * ($cluster_size{$spc}{$c});
		    }
		    print "Cluster $c\tSpecies $spc\t Size $cluster_size{$spc}{$c} \t In-density $indensity{$spc}{$c}\tOut-density $outdensity{$spc}{$c}\t Cluster modularity $total_cluster_modularity{$spc}{$c} \t Cluster modularity weighted $total_cluster_modularity_weighted{$spc}{$c} \t Cluster modularity by Size $total_cluster_modularity_bySize{$spc}{$c} \t Total in-cluster edge $total_cluster_edge{$spc}{$c}\t Total in-cluster degree $total_cluster_degree{$spc}{$c} \t Weighted in/out density ratio $weighted_in_out_ratio{$spc}{$c}\n";
		    $sumindensity += $cluster_size{$spc}{$c} * $indensity{$spc}{$c};
		    $countindensity += $cluster_size{$spc}{$c};
		    $sumoutdensity += $cluster_size{$spc}{$c} * $outdensity{$spc}{$c};
		    $countoutdensity += $cluster_size{$spc}{$c};
		    $sum_weighted_in_out_ratio += $weighted_in_out_ratio{$spc}{$c};
		}
    }
    $sumindensity /= $countindensity;
    $sumoutdensity /= $countoutdensity;
    my $out_in_ratio = $sumoutdensity/$sumindensity;
    # $sum_product_verDen_verIntro /= ($clusterNum * 3);
    print "All clusters on average have in-density of $sumindensity and out-density of $sumoutdensity and weighted ratio (out:in) ratio of and Sum of Modularity $sum_modularity_cluster_spc and Sum of Modularity Weighted $sum_modularity_cluster_spc_weighted and Sum of Modularity divided by Size $sum_modularity_cluster_spc_bySize and Weighted in/out density ratio $sum_weighted_in_out_ratio and Cluster num = $clusterNum\n";
}