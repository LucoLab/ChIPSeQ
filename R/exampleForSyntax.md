ChIP-SEQ Analysis
================
Jean-Philippe Villemin
May 22, 2017

Generate document body from comments
====================================

All the features from markdown and markdown supported within .Rmd documents, I was able to get from within R scripts. Here are some that I tested and use most frequently:

-   Smart comment fomatting in your R script generate the body and headers of the document
    -   Simply tweak your comments to begin with `#'` instead of just `#`
-   Create markdown headers as normal: `#' #` for h1, `#' ##` for h2, etc.
-   Add two spaces to the end of a comment line to start a new line (just like regular markdown)
-   Add `toc: true` to YAML frontmatter to create a table of contents with links like the one at the top of this page that links to h1, h2 & h3's indented like so:
    -   h1
        -   h2
            -   h3
-   Modify YAML to change syntax highlighting style (I'm using zenburn), author, title, theme, and all the good stuff you're used to setting in Rmd documents.
-   Sub-bullets like the ones above are created by a `#' *` with 4 spaces per level of indentation.
-   Surround text with `*` to *italicize*
-   Surround text with `**` to **bold**
-   Surround text with `***` to ***italicize & bold***
-   Skip lines with `#' <br>`
-   Keep comments in code, but hide from printing in report with `#' <!-- this text will not print in report -->`
-   Add hyperlinks:
    -   [Rmarkdown cheatsheet](http://rmarkdown.rstudio.com/RMarkdownCheatSheet.pdf)
    -   [Rmarkdown Reference Guide](http://rmarkdown.rstudio.com/RMarkdownReferenceGuide.pdf)
    -   [Compiling R notebooks from R Scripts](http://rmarkdown.rstudio.com/r_notebook_format.html)

The report begins here.

You can use the special syntax {{code}} to embed inline expressions, e.g. 2.1935703 is the mean of x plus 2. The code itself may contain braces, but these are not checked. Thus, perfectly valid (though very strange) R code such as `{{2 + 3}} - {{4 - 5}}` can lead to errors because `2 + 3}} - {{4 - 5` will be treated as inline code.

Now we continue writing the report. We can draw plots as well.

``` {.r}
plot(x)
```

![](figure/silk-test-b-1.png)

Actually you do not have to write chunk options, in which case knitr will use default options. For example, the code below has no options attached:

``` {.r}
var(x)
```

    ## [1] 0.6577564

``` {.r}
quantile(x)
```

    ##          0%         25%         50%         75%        100% 
    ## -0.56047565 -0.23017749  0.07050839  0.12928774  1.55870831

And you can also write two chunks successively like this:

``` {.r}
sum(x^2) # chi-square distribution with df 5
```

    ## [1] 2.818373

``` {.r}
sum((x - mean(x))^2) # df is 4 now
```

    ## [1] 2.631026

Done. Call spin('knitr-spin.R') to make silk from sow's ear now and knit a lovely purse.
