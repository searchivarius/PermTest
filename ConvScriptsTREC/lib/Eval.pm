#
# This code is released under the
# Apache License Version 2.0 http://www.apache.org/licenses/.
#
# (c) Leonid Boytsov, http://boytsov.info
#
sub Median {
    my $ArrRef = shift;

    my @a = sort { $a <=> $b } @$ArrRef;
    my $N = scalar(@a);
    my $m = int($N/2);
    if ($N % 2 != 0) {
        return $a[$m];
    }
    return undef if (!$N);

    return ($a[$m] + $a[$m - 1]) / 2;
}

1;
