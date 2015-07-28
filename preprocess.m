load('~/OneDrive/repos/safari-task-model_old/traindata.mat')

% tours table needs this
tours_subject = [];
tours_trial = [];
tours_sector = [];
tours_animal = [];
tours_question = [];
tours_tour = [];
tours_answer = [];
tours_optima = [];
tours_t_question = [];
tours_t_response = [];
tours_b_response = [];

% trials table needs this
trials_subject = [];
trials_trial = [];
trials_sector = [];
trials_length = [];
trials_animal = [];
trials_questions_sector = [];
trials_questions_biggersmaller = [];
trials_posteriors_new = [];
trials_answers_old = [];
trials_answers_new = [];
trials_posteriors_old = [];
trials_t_question = [];
trials_t_response = [];
trials_b_response = [];

prevsector = 0;
tour = 0;

for subject = 1:length(traindata)
	% disp 'Processing subject' , subject
	tours_sessions = find(traindata(subject).stimlist.tour_or_trials == 1);
	trials_sessions = find(traindata(subject).stimlist.tour_or_trials == 2);
	for session = tours_sessions
		% disp 'session' , session
		for sector = traindata(subject).stimlist.tours.sectors(session,:)
			% disp 'sector' , sector
			if prevsector ~= sector
				tour = tour + 1; 
			end
			prevsector = sector; 
			for trial = 1:length(traindata(subject).stimlist.tours.answers{session})
				% disp 'trial' , trial
				tours_subject = [tours_subject; subject ];
				tours_trial = [tours_trial; trial ];
				tours_sector = [tours_sector; sector ];
				tours_animal =	[tours_animal; traindata(subject).stimlist.tours.animals{session,sector}(trial)];
				tours_question = [tours_question; traindata(subject).stimlist.tours.questions{session,sector}(trial,:)];
				tours_tour = [tours_tour; tour];
				tours_answer = [tours_answer; traindata(subject).stimlist.tours.answers{session,sector}(trial)];
				tours_optima = [tours_optima; traindata(subject).stimlist.tours.optima{session,sector}(trial)];
				tours_t_question = [tours_t_question; traindata(subject).tours.t.question{session,sector}(trial)];
				tours_t_response = [tours_t_response; traindata(subject).tours.t.response{session,sector}(trial)];
				tours_b_response = [tours_b_response; traindata(subject).tours.b.response{session,sector}(trial)];
			end
		end
	end

	for session = trials_sessions
		% disp 'session' , session
		for sector = 1:length(traindata(subject).stimlist.trials.lengths{session})
			for trial = 1:traindata(subject).stimlist.trials.lengths{session}(sector)
				% disp 'trial' , trial
				trials_subject = [trials_subject ; subject];
				trials_trial = [trials_trial ; trial];
				trials_sector = [trials_sector ; sector];
				trials_length = [trials_length ; traindata(subject).stimlist.trials.lengths{session}(sector)];
				trials_animal = [trials_animal ; traindata(subject).stimlist.trials.animals{session}{sector}(trial)];
				trials_questions_sector = [trials_questions_sector ; traindata(subject).stimlist.trials.questions_sectors{session}(sector,:)];
				trials_questions_biggersmaller = [trials_questions_biggersmaller ; traindata(subject).stimlist.trials.questions_biggersmaller{session}(sector)];
				trials_posteriors_new = [trials_posteriors_new ; traindata(subject).stimlist.trials.posteriors_new{session}{sector}(trial,:)];
				trials_posteriors_old = [trials_posteriors_old ; traindata(subject).stimlist.trials.posteriors_old{session}{sector}(trial,:)];
				trials_answers_old = [trials_answers_old ; traindata(subject).stimlist.trials.answers_old{session}(sector)];
				trials_answers_new = [trials_answers_new ; traindata(subject).stimlist.trials.answers_new{session}(sector)];
				trials_t_question = [trials_t_question ; traindata(subject).trials.t.question{session}(trial)];
				trials_t_response = [trials_t_response ; traindata(subject).trials.t.response{session}(trial)];
				trials_b_response = [trials_b_response ; traindata(subject).trials.b.response{session}(trial)];
			end
		end
	end
end
% now create data tables
tours_varnames = {'subject', 'trial', 'sector', 'animal', 'question', 'tour', 'answer', 'optima', 't_question', 't_response', 'response'};
tours_tab = table(tours_subject, tours_trial, tours_sector, tours_animal, tours_question, tours_tour, tours_answer, tours_optima, tours_t_question, tours_t_response, tours_b_response, 'VariableNames', tours_varnames);
trials_varnames = {'subject', 'trial', 'sector', 'length', 'animal', 'questions_sector', 'questions_biggersmaller', 'posteriors_new', 'answers_old', 't_question', 't_response', 'response'}; 
trials_tab = table(trials_subject, trials_trial, trials_sector, trials_length, trials_animal, trials_questions_sector, trials_questions_biggersmaller, trials_posteriors_new, trials_answers_old, trials_t_question, trials_t_response, trials_b_response, 'VariableNames', trials_varnames);

writetable(tours_tab, 'tours.csv');
writetable(trials_tab, 'trials.csv');