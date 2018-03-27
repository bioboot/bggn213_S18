# Foundations of Bioinformatics course website  

Author and course director: Barry J. Grant (<bjgrant@ucsd.edu>)

Foundations of Bioinformatics (BGGN-213) is a UCSD bioscience graduate course run in the Spring quarter of 2018.  
For full detailes see: **https://bioboot.github.io/bggn213_S18/** 


**Overview:** This is a simple [jekell based static site](http://jekyllrb.com/docs/home/). To view locally on your own machine (i.e. before pushing or submitting a pull 
request to this [bioboot GitHub](https://github.com/bioboot/bimm143_S18) repo) 
you will need to have the **jekyll** and **github-pages** gem setup (see further 
below for full instructions)


### View, Edit and Contribute Changes Locally
This is a simple [jekell based static site](http://jekyllrb.com/docs/home/). 

To view locally on your own machine (i.e. before pushing or submitting a pull 
request to this [bioboot GitHub](https://github.com/bioboot/bimm143_W18/) repo) 
you will need to have the **jekyll** and **github-pages** gem setup, i.e.:

## Install ruby, jekyll and github-pages 
To be added - for now just use `homebrew` on mac.

Next consider updating RubyGems first (may need sudo for these).

	sudo gem update --system

Then install the Jekyll Gem and the GitHub Gem

	gem install jekyll
	gem install github-pages

Optional: Pygments python based syntax highlighter

	pip install Pygments


## Basics of Jekyll websites
Jekyll websites are configured based on the contents of the various underscore prefixed files and folders. You can find out more about these here: http://jekyllrb.com/docs/structure/

However, most likely you will want to leave most of these alone and just add  
content to the lectures.md file and create new files in the **class-material/** 
directory (i.e. add lecture slides, handouts, cheat-sheets etc.)

Please remember that all content is on the **gh-pages** branch! 
So you will want to be working on this branch and push back to this branch.

A typical workflow for folks that have been added as **"Collaborators"** would look something like this:

	## One time only clone
	git clone https://github.com/bioboot/bimm143_W18.git
	cd bimm143_W18

	## Edit your files (e.g. lectures.md, _data/authors.yml, _config.yml)
	vi lectures.md

	## Check changes localy
	jekyll serve

	## Pull recent changes
	git pull origin gh-pages

	## Stage, commit and push your changes
	git status
	git add lectures.md
	git commit -m "Your msg about changes"
	git push origin gh-pages



### Roll forward instructions...

To roll forward for a new years class follow the steps below (assuming you already have **jekyll** and **github-pages** 
setup on your local machine):

Git clone old site to a new dir

  	cd ~/Dropbox/Teaching
  	mkdir bggn213_2018
  	cd bggn213_2018
  	git clone git@github.com:bioboot/bggn213_f17.git bggn213_f18
  	cd bggn213_f18/
  
Update `_config.yml`, `_data/authors.yml`, `index.md` and of course `lectures.md`. 

In particular: 
- in `_config.yml` change all **bimm143_W18** to the new course name, including **title**.
- in `_data/authors.yml` change the **Course Syllabus** link, logo src, and **Instructional Assistant**.
- in  `index.md` update the bimm143_W18 etc links and names, syllabus link, initial **course evaluation** link 


Also remember to change the dates and the **pre-course questionnaire** and **post-course evaluation forms** as well as the lecture dates. Go through the regular `git add`, `git commit -m` cycle. But don’t yet push to GitHub (as we will want a new repo for this years class).
  

On GitHub make a new repo (Use the “+” sign and name it `bggn213_f18` to match your local directory name. This name matching is purely for convenience).

Then on the local machine change your remotes to point to this new repo.

  	git remotes -v   
  	git remote rm origin  

Now add our new repo and push changes:  

  	git remote add origin git@github.com:bioboot/bggn213_f18.git  
  	git push -u origin gh-pages  

Then preview your new site online: https://bioboot.github.io/web-2017/ and visit the repo itself to see if everything is ship-shape: https://github.com/bioboot/web-2017  
